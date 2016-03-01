from django.contrib.auth.decorators import login_required
import datetime
import csv
import json
import sys
from django.views.decorators.csrf import csrf_exempt
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404
from django.http import HttpResponse

from xbrowse.analysis_modules.combine_mendelian_families import get_variants_by_family_for_gene
from xbrowse_server.analysis.diagnostic_search import get_gene_diangostic_info
from xbrowse_server.base.models import Project, Family, FamilySearchFlag, VariantNote, ProjectTag, VariantTag
from xbrowse_server.api.utils import get_project_and_family_for_user, get_project_and_cohort_for_user, add_extra_info_to_variants_family
from xbrowse_server.api import utils as api_utils
from xbrowse_server.api import forms as api_forms
from xbrowse_server.mall import get_reference, get_datastore, get_mall
from xbrowse_server.search_cache import utils as cache_utils
from xbrowse_server.decorators import log_request
from xbrowse_server.server_utils import JSONResponse
from xbrowse.variant_search import cohort as cohort_search
from xbrowse import Variant
from xbrowse.analysis_modules.mendelian_variant_search import MendelianVariantSearchSpec
from xbrowse.core import displays as xbrowse_displays
from xbrowse_server import server_utils
from xbrowse_server import user_controls

import bigquery
import time

project_id = 'xbrowse-prod'

client = bigquery.get_client(project_id, json_key_file= "/Users/weisburd/bin/data-a1dbeee256bd.json", readonly=True)

#async query
#job_id, _results = client.query('SELECT * FROM dataset.my_table LIMIT 1000')
#complete, row_count = client.check_job(job_id)
# if complete:
#    results = client.get_query_rows(job_id)

#synchronous query
def query(q):
    print(q)
    start = time.clock()

    _job_id, results = client.query(q, timeout=10000)

    duration = (time.clock() - start)*1000

    return results, duration
    #except bigquery.BigQueryTimeoutException:



def compute_sql_query(project, family, search_spec):
    """{
        "gene_burden_filter": "",
        "genotype_inheritance_filter": "",
        "inheritance_mode": "x_linked_recessive",
        "quality_filter": {
            "min_ab": 25,
            "min_gq": 20,
            "vcf_filter": "pass"
        },
        "search_mode": "standard_inheritance",
        "variant_filter": {
            "genes": [
                "ENSG00000173991",
                "ENSG00000155657",
                "ENSG00000135636"
            ],
            "ref_freqs": [
                [
                    "1kg_wgs_phase3",
                    0.01
                ],
                [
                    "exac_v3",
                    0.01
                ]
            ]
        }
    }"""
    
    project_id_map = {'1kg_subset': '1kg'}
    
    project_id = project_id_map.get(project.project_id, project.project_id)
    
    #sql_query += " LEFT JOIN [%(project_id)s.variant_impacts_chosen_transcript] as vi ON v.xpos=vi.xpos AND v.ref=vi.ref AND v.alt=vi.alt " % locals()
    #sql_query += " LEFT JOIN [reference_data_for_variants.all_data_v2015_11_30] as rd ON v.xpos=rd.xpos AND v.ref=rd.ref AND v.alt=rd.alt " % locals()
    variant_genotypes_sql_query = """
SELECT
  xpos,
  chrom,
  pos,
  ref,
  alt,
  sample__NA19675__GT,
  sample__NA19678__GT,
  sample__NA19679__GT
FROM
  [%(project_id)s.variant_genotypes_pivoted]
""" % locals()

    variant_genotypes_where_conditions_list = []
    
    search_spec = search_spec.toJSON()
    print(json.dumps(search_spec))
    
    inheritance_mode = search_spec['inheritance_mode']
    
    num_alt = {}
    genotype_filter = {}
    
    if inheritance_mode == "recessive":  # TODO fix this and compound het search
        inheritance_mode = 'homozygous_recessive'   
        
    if inheritance_mode in ('homozygous_recessive', 'x_linked_recessive'):
        
        if inheritance_mode == 'x_linked_recessive':
            variant_genotypes_where_conditions_list.append("chrom='X'")
        
        # set affected and unaffected genotypes first, without inheritance
        for indiv_id, individual in family.individuals.items():
            print("#### " + indiv_id)
            if individual.affected_status == 'affected':
                genotype_filter[indiv_id] = "sample__%(indiv_id)s__num_alt = 2" % locals()   #'alt_alt'
            elif individual.affected_status == 'unaffected':
                genotype_filter[indiv_id] = "sample__%(indiv_id)s__num_alt <= 1" % locals()  #'has_ref'           
    
        
        # now account for parental relationships (but no others)
        for indiv_id, individual in family.individuals.items():    
            
            if individual.affected_status != 'affected':
                continue
            
            father = family.get_individual(individual.paternal_id)
            mother = family.get_individual(individual.maternal_id)

            if mother is not None:
                if mother.affected_status == 'unaffected':
                    genotype_filter[individual.maternal_id] = "sample__%s__num_alt = 1" % individual.maternal_id  #'ref_alt'

            if father is not None:
                # if father is unaffected, should be ref_ref instead of ref_alt
                if father.affected_status == 'unaffected':
                    if inheritance_mode == 'homozygous_recessive':
                        num_alt = 1 #'ref_alt'
                    elif inheritance_mode == 'x_linked_recessive':
                        num_alt = 0 #'ref_ref'    
                    else:
                        raise Exception("Unexpected recessive inheritance mode: %s" % inheritance_mode)
                    
                    genotype_filter[individual.paternal_id] = "sample__%s__num_alt = %d" % (individual.paternal_id, num_alt)
                     
    elif inheritance_mode in ('de_novo', 'dominant'):
        # set affected and unaffected genotypes first, without inheritance
        for indiv_id, individual in family.individuals.items():
            if individual.affected_status == 'affected':
                if inheritance_mode == 'dominant':
                    genotype_filter[indiv_id] = "sample__%(indiv_id)s__num_alt >= 1" % locals()  #'ref_alt'
                elif inheritance_mode == 'de_novo':
                    genotype_filter[indiv_id] = "sample__%(indiv_id)s__num_alt = 1" % locals()  #'ref_alt'
                else:
                    raise Exception("Unexpected dominant inheritance mode: %s" % inheritance_mode)
            elif individual.affected_status == 'unaffected':
                genotype_filter[individual.paternal_id] = "sample__%(indiv_id)s__num_alt = 0" % locals()  #'ref_ref'

    variant_genotypes_where_conditions_list.extend(genotype_filter.values())
    
    # add quality filters
    quality_filter = {}
    if 'quality_filter' in search_spec:
        quality_filter = search_spec['quality_filter']
        min_ab = int(quality_filter['min_ab'])/100.
        min_gq = int(quality_filter['min_gq'])
        #vcf_filter = quality_filter['vcf_filter']
        for indiv_id, individual in family.individuals.items():
            variant_genotypes_where_conditions_list.append( "(sample__%(indiv_id)s__AB is null or sample__%(indiv_id)s__AB > %(min_ab)s)" % locals() ) #'ref_ref'
            variant_genotypes_where_conditions_list.append( "sample__%(indiv_id)s__GQ > %(min_gq)s" % locals() )  #'ref_ref'
    
    if variant_genotypes_where_conditions_list:
        variant_genotypes_sql_query += "WHERE " + " and ".join(variant_genotypes_where_conditions_list)    
    
    variant_genotypes_where_conditions_list = None
    
    
    variant_impacts_where_conditions_list = []
    variant_impacts_sql_query = """
SELECT
    xpos,
    ref,
    alt,
    gene,
    feature,
    consequence,
    severity_rank,
    hgvsc,
    hgvsp,
    symbol,
    symbol_source,
    biotype,
    canonical,
    lof_info,
    lof_flags,
    lof_filter,
    lof,
    sift_pred,
    polyphen2_hvar_pred,
    cadd_phred,
    mutationtaster_pred,
    fathmm_pred,
    metasvm_pred
  FROM
    [%(project_id)s.variant_impacts_chosen_transcript]
    """ % locals()

    if 'variant_filter' in search_spec: 
        variant_filter = search_spec['variant_filter']
        #if 'ref_freqs' in variant_filter:
        #    for ref_freq in variant_filter['ref_freqs']:
        #        # ["1kg_wgs_phase3", 0.01], ["1kg_wgs_phase3_popmax", 0.01], ["exac_v3", 0.01], ["exac_v3_popmax", 0.01]
        #        if ref_freq[0] == '1kg_wgs_phase3':
        #            variant_impacts_where_conditions_list.append('(rd.g1k_AF is null or rd.g1k_AF <= %s)' % ref_freq[1])
        #        elif ref_freq[0] == '1kg_wgs_phase3_popmax':
        #            variant_impacts_where_conditions_list.append('(rd.g1k_AF_POPMAX is null or rd.g1k_AF_POPMAX <= %s)' % ref_freq[1])
        #        elif ref_freq[0] == 'exac_v3':
        #            variant_impacts_where_conditions_list.append('(rd.exac_v3_AF is null or rd.exac_v3_AF <= %s)' % ref_freq[1])
        #        elif ref_freq[0] == 'exac_v3_popmax':
        #            variant_impacts_where_conditions_list.append('(float(rd.exac_v3_AF_POPMAX) is null or float(rd.exac_v3_AF_POPMAX) <= %s)' % ref_freq[1])
        if 'genes' in variant_filter:
            variant_impacts_where_conditions_list.append("gene in ('%s')" % str("','".join([gene for gene in variant_filter['genes']])))
            
        if 'so_annotations' in variant_filter:
            variant_impacts_where_conditions_list.append("consequence in ('%s')" % str("','".join([so_annot for so_annot in variant_filter['so_annotations']])))
        
    if variant_impacts_where_conditions_list:
        variant_impacts_sql_query += "WHERE " + " and ".join(variant_impacts_where_conditions_list)    
    
    sql_query = """
SELECT
  v.chrom as chrom,
  v.pos as pos,
  v.ref as ref,
  v.alt as alt, 
  sample__NA19675__GT,
  sample__NA19678__GT,
  sample__NA19679__GT,
  gene,
  consequence,
  hgvsc,
  hgvsp,
  symbol,
  symbol_source,
  biotype,
  canonical,
  lof_info,
  lof_flags,
  lof_filter,
  lof,
  sift_pred,
  polyphen2_hvar_pred,
  cadd_phred,
  mutationtaster_pred,
  fathmm_pred,
  metasvm_pred
FROM (%(variant_genotypes_sql_query)s) as v 
JOIN (%(variant_impacts_sql_query)s) as vi
ON v.xpos=vi.xpos AND v.ref=vi.ref AND v.alt=vi.alt
""" % locals()
    return sql_query
        
        

@login_required
#@log_request('variant_search_api')
def variant_search_api(request):
    #sys.stderr.write(__file__ + " #### variant_search_api(" + str(request) + ")\n")
    
    project, family = get_project_and_family_for_user(request.user, request.GET)
    
    form = api_forms.MendelianVariantSearchForm(request.GET)
    if form.is_valid():

        search_spec = form.cleaned_data['search_spec']
        sql_query = compute_sql_query(project, family.xfamily(), search_spec)
        print("Running sql query: \n############\n" + sql_query + "\n\n############")
        result, duration = query(sql_query)
        print("Duration: %s millisec" % str(duration))
        print_header = True
        for r in result:
            if print_header:
                print("\t".join(k for k in sorted(r.keys())))
                print("\t".join('-'*len(k) for k in sorted(r.keys())))
                print_header = False
            print("\t".join(str(r[k]) for k in sorted(r.keys())))
                
            
        
        search_spec.family_id = family.family_id
        
        variants = api_utils.calculate_mendelian_variant_search(search_spec, family.xfamily())
        search_hash = cache_utils.save_results_for_spec(project.project_id, search_spec.toJSON(), [v.toJSON() for v in variants])
        add_extra_info_to_variants_family(get_reference(), family, variants)

        return_type = request.GET.get('return_type', 'json')
        if return_type == 'json':
            return JSONResponse({
                'is_error': False,
                'variants': [v.toJSON() for v in variants],
                'search_hash': search_hash,
            })
        elif return_type == 'csv':
            return ''
        else:
            return HttpResponse("Return type not implemented")

    else:
        return JSONResponse({
            'is_error': True,
            'error': server_utils.form_error_string(form)
        })


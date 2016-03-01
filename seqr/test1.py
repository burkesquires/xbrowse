"""


Service account key
Created new service account: xbrowse@xbrowse-prod.iam.gserviceaccount.com
Downloaded json key: ./bin/data-a1dbeee256bd.json

https://github.com/tylertreat/BigQuery-Python
"""

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


result, duration = query("""SELECT
  i.vep.symbol,
  v.chrom,
  v.pos,
  v.ref,
  v.alt,
  v.sample__100AO_HS_1__num_alt
FROM
  [INMR_v9.variant_genotypes_pivoted] AS v
JOIN EACH [INMR_v9.variant_impacts] AS i
ON
  v.xpos=i.xpos
  AND v.ref=i.ref
  AND v.alt=i.alt
WHERE
v.sample__100AO_HS_1__num_alt=2 AND
v.sample__102AP_BG_1__num_alt=1 AND
v.sample__103AP_PG_1__num_alt=1
""")

print("query returned %s results and took %s millisec" % (len(result), duration))
print(result[0])




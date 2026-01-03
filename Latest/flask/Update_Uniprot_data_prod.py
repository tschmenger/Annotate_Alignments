import psycopg2
import json
import os
from psycopg2 import sql
WURZELDIR = ""
#########################################################################################################################################################################################################################
## turn the following line on when on production
os.environ.setdefault("URL_PREFIX", "/CALVI")

## turn the following two lines off when on production
#os.environ.setdefault("URL_PREFIX", "")
#os.environ.setdefault("LOCAL_MODE", "true")

###
URL_PREFIX = os.getenv("URL_PREFIX", "")  # Will be "" locally, and "/CALVI" in production
                                          # export URL_PREFIX="/CALVI"
LOCAL_MODE = os.getenv("LOCAL_MODE", "false").lower() == "true"
                                          # export LOCAL_MODE=true
###

if LOCAL_MODE:
    WURZELDIR = "/mnt/c/Users/TSchm/Desktop/Coding/Projekt 3 - Alignments"
else:
    sourcefilepath = "/home/calvi/"
#########################################################################################################################################################################################################################
def fetch_uniprot_rows(cur, primacc_list):
    if not primacc_list:
        return {}

    query = sql.SQL("""
        SELECT primacc, genus, modus, position
        FROM uniprotinfo
        WHERE primacc IN %s
    """)

    cur.execute(query, (tuple(primacc_list),))
    rows = cur.fetchall()

    return {
        primacc.strip(): (genus, modus, position)
        for primacc, genus, modus, position in rows
    }

test_batch = ["P61586", "Q12888 "]
test_batch = [x.strip() for x in test_batch]

with open('../config.json') as f:
    config = json.load(f)

conn = psycopg2.connect(
    dbname=config["DB_NAME"],
    user="postgres"
)
cur = conn.cursor()

before = fetch_uniprot_rows(cur, test_batch)



# Create staging table
cur.execute("""
CREATE TEMP TABLE IF NOT EXISTS uniprotinfo_stage (
    primacc  TEXT,
    genus    TEXT,
    modus    TEXT,
    position TEXT
);
""")

# Load file
if LOCAL_MODE:
    file_path = os.path.join(WURZELDIR, "psql_DB", "UniprotData_V2_20260102.txt")
else:
    file_path = os.path.join(sourcefilepath, "psql_DB", "UniprotData_V2_20260102.txt")
  
with open(file_path, "r") as f:
    cur.copy_from(
        f,
        "uniprotinfo_stage",
        sep="\t",
        columns=("primacc", "genus", "modus", "position")
    )

# Trim whitespace
cur.execute("""
UPDATE uniprotinfo_stage
SET primacc = TRIM(primacc),
    genus = TRIM(genus),
    modus = TRIM(modus),
    position = TRIM(position);
""")


# UPSERT
cur.execute("""
DELETE FROM uniprotinfo t
WHERE NOT EXISTS (
    SELECT 1
    FROM uniprotinfo_stage s
    WHERE s.primacc = t.primacc
      AND s.genus = t.genus
      AND s.modus = t.modus
      AND s.position = t.position
);
""")

cur.execute("""
INSERT INTO uniprotinfo (primacc, genus, modus, position)
SELECT s.primacc, s.genus, s.modus, s.position
FROM uniprotinfo_stage s
LEFT JOIN uniprotinfo t
  ON t.primacc = s.primacc
 AND t.genus = s.genus
 AND t.modus = s.modus
 AND t.position = s.position
WHERE t.primacc IS NULL;
""")


conn.commit()

after = fetch_uniprot_rows(cur, test_batch)

cur.close()
conn.close()



print(before)
print(after)


from pathlib import Path
import csv, urllib.request

manifest=Path("data/public_data_manifest.csv")
out=Path("data/public"); out.mkdir(parents=True,exist_ok=True)
for row in csv.DictReader(manifest.open()):
    suffix=".csv"
    target=out/(row["name"]+suffix)
    print("fetch",row["name"],"->",target)
    urllib.request.urlretrieve(row["url"],target)

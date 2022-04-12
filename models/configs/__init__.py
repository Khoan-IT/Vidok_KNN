import os
import json
from collections import defaultdict

current_dir = os.path.dirname(__file__)
config_files = os.listdir(current_dir)
config_files = list(filter(lambda x: x.endswith(".json"), config_files))

model_configs = defaultdict(lambda x: None)

for cfile in config_files:
    model_name = cfile.split('.')[0]
    model_configs[model_name] = json.load(
                                open(os.path.join(current_dir, cfile), "r", encoding="utf-8")
                                        )

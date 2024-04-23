# Copyright 2024 Huawei Technologies Co., Ltd.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# author: Pal Andras Papp, Georg Anegg

import os
from dataclasses import dataclass
from sys import platform

per = "/" if platform == "linux" else "\\"
merge_dir = f"..{per}temp_merged_instances"

@dataclass
class Settings:
    NUMA: bool = True
    Contract: bool = False
    TimeLimit: int = 3

@dataclass
class Config:
    file_mode: str = "Dir"
    location: str = f"..{per}instances"
    settings: Settings = Settings()


def read_config(filename):
    config = Config()
    with open(filename, 'r') as file:
        for line in file:
            cleaned = line.strip().split()

            if line.startswith('%') or len(cleaned) == 0:
                continue

            if len(cleaned) >= 2 and cleaned[1] == "=":
                cleaned.pop(1)
            if len(cleaned) < 2:
                print("Config file format error. Applying default configs instead.")
                return Config()

            if cleaned[0] in ["File", "Dir", "Combine"]:
                config.file_mode = cleaned[0]
                if cleaned[0] == "Combine":
                    if len(cleaned) >= 3:
                        config.location = [cleaned[1], cleaned[2]]
                        if not os.path.isdir(cleaned[1]) or not os.path.isdir(cleaned[2]):
                            print("Config error (incorrect file/directory). Applying default configs instead.")
                            return Config()
                    else:
                        print("Config file format error. Applying default configs instead.")
                        return Config()
                else:
                    config.location = cleaned[1]
                    if (cleaned[0] == "File" and not os.path.isfile(cleaned[1])) or (
                            cleaned[0] == "Dir" and not os.path.isdir(cleaned[1])):
                        print("Config error (incorrect file/directory). Applying default configs instead.")
                        return Config()

            if cleaned[0] in ["NUMA", "Contract"]:
                if cleaned[1] not in ["True", "False", "Yes", "No", "1", "0"]:
                    print("Config file format error. Applying default configs instead.")
                    return Config()
                value = True if cleaned[1] in ["True", "Yes", "1"] else False
                if cleaned[0] == "NUMA":
                    config.settings.NUMA = value
                else:
                    config.settings.Contract = value

            if cleaned[0] == "TimeLimit":
                limit = int(cleaned[1])
                if limit < 1:
                    print("Config file format error. Applying default configs instead.")
                    return Config()
                config.settings.TimeLimit = limit

    return config


def CopyFile(oldname, newname):
    if not os.path.exists(oldname):
        print(f'Unable to copy {oldname}: not found.')
        return False
    oldfile = open(oldname, 'r')
    newfile = open(newname, 'w')
    newfile.write(oldfile.read())
    oldfile.close()
    newfile.close()
    return True


def GetFileNames(settings):
    if settings.file_mode == "File":
        names = [settings.location]
    elif settings.file_mode == "Dir":
        names = [f'{settings.location}{per}{file}' for file in os.listdir(settings.location)]
    else:
        names = [f'{merge_dir}{per}{file}' for file in os.listdir(merge_dir)]
    for name in names:
        if not os.path.isfile(name):
            names.remove(name)
    return names


def CreateMergedFiles(config):
    dir_DAGs = os.listdir(config.location[0])
    dir_Params = os.listdir(config.location[1])
    if not os.path.isdir(merge_dir):
        os.mkdir(merge_dir)

    for file_DAG in dir_DAGs:
        for file_Params in dir_Params:
            if file_DAG.endswith(".txt"):
                new_name = file_DAG[:-4] + '_' + file_Params
            else:
                new_name = file_DAG + file_Params
            CopyFile(f'{config.location[0]}{per}{file_DAG}', f'{merge_dir}{per}{new_name}')
            base_file = open(f'{merge_dir}{per}{new_name}', 'a+')
            to_append = open(f'{config.location[1]}{per}{file_Params}', 'r')
            base_file.write(to_append.read())
            base_file.close()
            to_append.close()


def RemoveMergedFiles():
    for file in os.scandir(merge_dir):
        if file.is_file():
            os.remove(f'{merge_dir}{per}{file.name}')

    if os.path.exists(merge_dir):
        os.rmdir(merge_dir)




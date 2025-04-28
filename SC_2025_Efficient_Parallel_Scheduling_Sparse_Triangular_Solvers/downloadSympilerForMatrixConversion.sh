#!/bin/bash

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

# @author Toni Boehnlein, Pal Andras Papp, Raphael S. Steiner

echo "This script will clone the sympiler/aggregation repository."
echo "Please take note that sympiler/aggregation is under the MIT license."
echo " "
echo "Please ensure the download you initiate is in line with applicable terms of use, "
echo "laws, and regulations."
echo " "
echo "On confirmation, the repository will be cloned and a patch will be applied to it."
echo "The patch will modify the code to allow for producing the permuted matrices of the suite sparse data set."
echo ""
echo "The script will also attempt to build the repository which requires cmake, make and METIS."
echo "For more information on the requirements, please refer to the README file in the cloned repository."
echo ""
read -p "I have taken note and agree [yes/no] " -r
echo ""
if [[ "$REPLY" = "yes" ]]; then
    echo "Cloning sympiler/aggregation repository..."
    git clone https://github.com/sympiler/aggregation.git

    patch -p0 < patches/permute.patch

    cd aggregation
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j$(nproc)

else
    echo "Exiting script."
    exit 1
fi
exit 0


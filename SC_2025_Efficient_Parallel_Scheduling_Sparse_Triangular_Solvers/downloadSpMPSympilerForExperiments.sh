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

echo "This script will clone the SpMP repository."
echo "Please take note of the SpMP license."
echo " "
echo "Copyright (c) 2015, Intel Corporation. All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Intel Corporation nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL INTEL CORPORATION BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
echo " "
echo "This script will also clone the sympiler/aggregation repository."
echo "Please take note that sympiler/aggregation is under the MIT license."
echo " "
echo "Please ensure the download you initiate is in line with applicable terms of use, "
echo "laws, and regulations."
echo " "
echo "On confirmation, the repositories will be cloned and a patch will be applied to the sympiler code."
echo "The patch will modify the code to allow for running SpTrSV experiments."
echo ""
echo "The script will also attempt to build SpMP and sympiler which requires cmake, make, and gcc"
echo "For more information on the requirements, please refer to the README file in the cloned repository."
echo ""
read -p "I have taken note and agree [yes/no] " -r
echo ""
if [[ "$REPLY" = "yes" ]]; then

    echo "Cloning IntelLabs/SpMP repository..."
    git clone https://github.com/IntelLabs/SpMP.git

    cd SpMP
    make CUSTOM_CXX=gcc -j$(nproc)
    export SPMPROOT=$(pwd)
    cd ..

    rm -rf aggregation

    echo "Cloning sympiler/aggregation repository..."
    git clone https://github.com/sympiler/aggregation.git

    patch -p0 < patches/experiments.patch

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


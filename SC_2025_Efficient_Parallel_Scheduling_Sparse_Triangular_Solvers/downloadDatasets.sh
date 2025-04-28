#!/bin/bash

#
#   Copyright 2021 Huawei Technologies Co., Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# Notice:
# This file has been modified from the original (https://github.com/Algebraic-Programming/ALP/blob/develop/tools/downloadDatasets.sh)
# in order to download the dataset from SuiteSparse used in a SC25 submission
# Modified by Toni Boehnlein and Raphael S. Steiner, date 24. April 25

downloadSS () {
	if [ ! -f ${1}.mtx ]; then
		if [ ! -f ${1}.tar.gz ]; then
			wget ${2} || exit 1
		fi
		if [ ! -f ${1}/${1}.mtx ]; then
			tar xf ${1}.tar.gz || exit 1
			mv ${1}/${1}.mtx ${1}.mtx || exit 1
			rm ${1}.tar.gz || exit 1
			rm -r ${1}/ || exit 1
		fi
	fi
}

OSP_DIR="$(pwd)/one-stop-parallel"

DATASETS_DIR="$(pwd)/data_set"
ALLMATRICES_DIR="$(pwd)/data_set/matrices"

SUITESPARSE_DIR="$(pwd)/data_set/suite_sparse"
METIS_DIR="$(pwd)/data_set/metis"
ICHOL_DIR="$(pwd)/data_set/ichol"
ERDOS_RENYI_DIR="$(pwd)/data_set/erdos_renyi"
NARROW_BANDWIDTH_DIR="$(pwd)/data_set/narrow_bandwidth"

SUITESPARSE_SET="af_0_k101 af_shell7 apache2 audikw_1 bmw7st_1 bmwcra_1 bone010 boneS01 boneS10 Bump_2911 bundle_adj consph Dubcova3 ecology2 Emilia_923 Fault_639 Flan_1565 G3_circuit Geo_1438 hood Hook_1498 inline_1 ldoor msdoor offshore parabolic_fem PFlow_742 Queen_4147 s3dkt3m2 Serena shipsec1 StocF-1465 thermal2"
ICHOL_SET="af_0_k101 af_shell7 apache2 audikw_1 bmw7st_1 bmwcra_1 bone010 boneS01 boneS10 Bump_2911 consph Dubcova3 ecology2 Emilia_923 Fault_639 Flan_1565 G3_circuit Geo_1438 hood Hook_1498 inline_1 ldoor msdoor offshore parabolic_fem PFlow_742 Queen_4147 s3dkt3m2 Serena shipsec1 StocF-1465 thermal2"
METIS_SET="af_0_k101 af_shell10 apache2 audikw_1 bmwcra_1 bone010 boneS10 bundle_adj cant consph crankseg_2 ecology2 Emilia_923 Fault_639 Flan_1565 G3_circuit Geo_1438 gyro hood Hook_1498 inline_1 ldoor m_t1 msdoor nasasrb PFlow_742 pwtk raefsky4 ship_003 shipsec8 StocF-1465 thermal2 tmt_sym x104"

SYMPILER_METIS_BIN="$(pwd)/aggregation/build/example/Hdagg_SpTRSV"

echo "This script will download matrices from the SuiteSparse Matrix Collection [1], which "
echo "is maintained by Tim Davis, Yifan Hu, and Scott Kolodziej."
echo " "
echo "It will further take these matrices and perform an incomplete Cholesky from Eigen and a Node permute from Metis on them."
echo " "
echo "Furthermore, 60 matrices will be sampled from six matrix distributions."
echo " "
echo "The matrices downloaded from the SuiteSparse Matrix Collection are:"
echo " - af_0_k101, Structural Problem"
echo " - af_shell7, Subsequent Structural Problem"
echo " - apache2, Structural Problem"
echo " - audikw_1, Structural Problem"
echo " - bmw7st_1, Structural Problem"
echo " - bmwcra_1, Structural Problem"
echo " - bone010, Model Reduction Problem"
echo " - boneS01, Model Reduction Problem"
echo " - boneS10, Model Reduction Problem"
echo " - Bump_2911, 2D/3D Problem [2, 3]"
echo " - bundle_adj, Computer Vision Problem"
echo " - consph, 2D/3D Problem [4]"
echo " - Dubcova3, 2D/3D Problem"
echo " - ecology2, 2D/3D Problem [5]"
echo " - Emilia_923, Structural Problem [6, 7]"
echo " - Fault_639, Structural Problem [8, 9 ,10, 11]"
echo " - Flan_1565, Structural Problem [11, 12]"
echo " - G3_circuit, Circuit Simulation Problem"
echo " - Geo_1438, Structural Problem [11]"
echo " - hood, Structural Problem"
echo " - Hook_1498, Structural Problem"
echo " - inline_1, Structural Problem"
echo " - ldoor, Structural Problem"
echo " - msdoor, Structural Problem [13]"
echo " - offshore, Electromagnetics Problem [14]"
echo " - parabolic_fem, Computational Fluid Dynamics Problem"
echo " - PFlow_742, 2D/3D Problem [3]"
echo " - Queen_4147, 2D/3D Problem [3]"
echo " - s3dkt3m2, Structural Problem"
echo " - Serena, Structural Problem [6]"
echo " - shipsec1, Structural Problem"
echo " - StocF-1465, Computational Fluid Dynamics Problem [7, 15]"
echo " - thermal2, Thermal Problem"
echo " - af_shell10, Structural Problem"
echo " - cant, 2D/3D Problem [4]"
echo " - crankseg_2, Structural Problem"
echo " - gyro, Model Reduction Problem"
echo " - m_t1, Structural Problem"
echo " - nasasrb, Structural Problem"
echo " - pwtk, Structural Problem"
echo " - raefsky4, Structural Problem"
echo " - ship_003, Structural Problem"
echo " - shipsec8, Structural Problem"
echo " - tmt_sym, Electromagnetics Problem"
echo " - x104, Structural Problem"
echo " "
echo "[1]  Timothy A. Davis and Yifan Hu. 2011. The University of Florida Sparse Matrix"
echo "     Collection."
echo "     ACM Transactions on Mathematical Software 38, 1, Article 1."
echo "[2]  C. Janna, M. Ferronato, G. Gambolati. \"Enhanced Block FSAI preconditioning using Domain Decomposition techniques\"."
echo "     SIAM Journal on Scientific Computing, 35, pp. S229-S249, 2013."
echo "[3]  C. Janna, M. Ferronato, G. Gambolati. \"The use of supernodes in factored sparse approximate inverse preconditioning\"."
echo "     SIAM Journal on Scientific Computing, submitted."
echo "[4]  S. Williams, L. Oliker, R. Vuduc, J. Shalf, K. Yelick, J. Demmel,"
echo "     \"Optimization of Sparse Matrix-Vector Multiplication on Emerging Multicore Platforms\","
echo "     Parallel Computing Volume 35, Issue 3, March 2009, Pages 178-194."
echo "     Special issue on Revolutionary Technologies for Acceleration of Emerging Petascale Applications."
echo "[5]  Brad McRae, National Center for Ecological Analysis and Synthesis Santa Barbara, CA."
echo "[6]  M. Ferronato, G. Gambolati, C. Janna, P. Teatini."
echo "     \"Geomechanical issues of anthropogenic CO2 sequestration in exploited gas fields\","
echo "     Energy Conversion and Management, 51, pp. 1918-1928, 2010."
echo "[7]  C. Janna, M. Ferronato."
echo "     \"Adaptive pattern research for block FSAI preconditionig\"."
echo "     SIAM Journal on Scientific Computing, to appear."
echo "[8]  M. Ferronato, G. Gambolati, C. Janna, P. Teatini."
echo "     \"Numerical modelling of regional faults in land subsidence prediction above gas/oil reservoirs\","
echo "     International Journal for Numerical and Analytical Methods in Geomechanics, 32, pp. 633-657, 2008."
echo "[9]  M. Ferronato, C. Janna, G. Gambolati."
echo "     \"Mixed constraint preconditioning in computational contact mechanics\","
echo "     Computer Methods in Applied Mechanics and Engineering, 197, pp. 3922-3931, 2008."
echo "[10] C. Janna, M. Ferronato, G. Gambolati."
echo "     \"Multilevel incomplete factorizations for the iterative solution of non-linear FE problems\"."
echo "     International Journal for Numerical Methods in Engineering, 80, pp. 651-670, 2009."
echo "[11] C. Janna, M. Ferronato, G. Gambolati."
echo "     \"A Block FSAI-ILU parallel preconditioner for symmetric positive definite linear systems\"."
echo "     SIAM Journal on Scientific Computing, 32, pp. 2468-2484, 2010."
echo "[12] C. Janna, A. Comerlati, G. Gambolati."
echo "     \"A comparison of projective and direct solvers for finite elements in elastostatics\"."
echo "     Advances in Engineering Software, 40, pp. 675-685, 2009."
echo "[13] from Parasol, http://www.parallab.uib.no/projects/parasol/data"
echo "[14] Evan Um, Geophysics, Stanford."
echo "[15] M. Ferronato, C. Janna, G. Pini."
echo "     \"Shifted FSAI preconditioners for the efficient parallel solution of non-linear groundwater flow models\"."
echo "     International Journal for Numerical Methods in Engineering, to appear."
echo " "
echo "Please take note of the attributions to SuiteSparse."
echo " "
echo "Please ensure the download you initiate is in line with applicable terms of use, "
echo "laws, and regulations."
echo " "
echo "Please use this script once and keep the datasets for future use."
echo ""
echo "On confirmation, the datasets will be downloaded to: ${DATASETS_DIR}"
echo ""
read -p "I have taken note and agree [yes/no] " -r
echo ""
if [[ "$REPLY" = "yes" ]]; then
	if [[ ! -d "${DATASETS_DIR}" ]]; then
		mkdir "${DATASETS_DIR}" || exit 1
	fi
	cd "${DATASETS_DIR}" || exit 1
	if [[ ! -d "${ALLMATRICES_DIR}" ]]; then
		mkdir "${ALLMATRICES_DIR}" || exit 1
	fi
	cd "${ALLMATRICES_DIR}" || exit 1
	downloadSS "af_0_k101" "https://suitesparse-collection-website.herokuapp.com/MM/Schenk_AFE/af_0_k101.tar.gz"
	downloadSS "af_shell7" "https://suitesparse-collection-website.herokuapp.com/MM/Schenk_AFE/af_shell7.tar.gz"
	downloadSS "apache2" "https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/apache2.tar.gz"
	downloadSS "audikw_1" "https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/audikw_1.tar.gz"
	downloadSS "bmw7st_1" "https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/bmw7st_1.tar.gz"
	downloadSS "bmwcra_1" "https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/bmwcra_1.tar.gz"
	downloadSS "bone010" "https://suitesparse-collection-website.herokuapp.com/MM/Oberwolfach/bone010.tar.gz"
	downloadSS "boneS01" "https://suitesparse-collection-website.herokuapp.com/MM/Oberwolfach/boneS01.tar.gz"
	downloadSS "boneS10" "https://suitesparse-collection-website.herokuapp.com/MM/Oberwolfach/boneS10.tar.gz"
	downloadSS "Bump_2911" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/Bump_2911.tar.gz"
	downloadSS "bundle_adj" "https://suitesparse-collection-website.herokuapp.com/MM/Mazaheri/bundle_adj.tar.gz"
	downloadSS "consph" "https://suitesparse-collection-website.herokuapp.com/MM/Williams/consph.tar.gz"
	downloadSS "Dubcova3" "https://suitesparse-collection-website.herokuapp.com/MM/UTEP/Dubcova3.tar.gz"
	downloadSS "ecology2" "https://suitesparse-collection-website.herokuapp.com/MM/McRae/ecology2.tar.gz"
	downloadSS "Emilia_923" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/Emilia_923.tar.gz"
	downloadSS "Fault_639" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/Fault_639.tar.gz"
	downloadSS "Flan_1565" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/Flan_1565.tar.gz"
	downloadSS "G3_circuit" "https://suitesparse-collection-website.herokuapp.com/MM/AMD/G3_circuit.tar.gz"
	downloadSS "Geo_1438" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/Geo_1438.tar.gz"
	downloadSS "hood" "https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/hood.tar.gz"
	downloadSS "Hook_1498" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/Hook_1498.tar.gz"
	downloadSS "inline_1" "https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/inline_1.tar.gz"
	downloadSS "ldoor" "https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/ldoor.tar.gz"
	downloadSS "msdoor" "https://suitesparse-collection-website.herokuapp.com/MM/INPRO/msdoor.tar.gz"
	downloadSS "offshore" "https://suitesparse-collection-website.herokuapp.com/MM/Um/offshore.tar.gz"
	downloadSS "parabolic_fem" "https://suitesparse-collection-website.herokuapp.com/MM/Wissgott/parabolic_fem.tar.gz"
	downloadSS "PFlow_742" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/PFlow_742.tar.gz"
	downloadSS "Queen_4147" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/Queen_4147.tar.gz"
	downloadSS "s3dkt3m2" "https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/s3dkt3m2.tar.gz"
	downloadSS "Serena" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/Serena.tar.gz"
	downloadSS "shipsec1" "https://suitesparse-collection-website.herokuapp.com/MM/DNVS/shipsec1.tar.gz"
	downloadSS "StocF-1465" "https://suitesparse-collection-website.herokuapp.com/MM/Janna/StocF-1465.tar.gz"
	downloadSS "thermal2" "https://suitesparse-collection-website.herokuapp.com/MM/Schmid/thermal2.tar.gz"
	downloadSS "af_shell10" "https://suitesparse-collection-website.herokuapp.com/MM/Schenk_AFE/af_shell10.tar.gz"
	downloadSS "cant" "https://suitesparse-collection-website.herokuapp.com/MM/Williams/cant.tar.gz"
	downloadSS "crankseg_2" "https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/crankseg_2.tar.gz"
	downloadSS "gyro" "https://suitesparse-collection-website.herokuapp.com/MM/Oberwolfach/gyro.tar.gz"
	downloadSS "m_t1" "https://suitesparse-collection-website.herokuapp.com/MM/DNVS/m_t1.tar.gz"
	downloadSS "nasasrb" "https://suitesparse-collection-website.herokuapp.com/MM/Nasa/nasasrb.tar.gz"
	downloadSS "pwtk" "https://suitesparse-collection-website.herokuapp.com/MM/Boeing/pwtk.tar.gz"
	downloadSS "raefsky4" "https://suitesparse-collection-website.herokuapp.com/MM/Simon/raefsky4.tar.gz"
	downloadSS "ship_003" "https://suitesparse-collection-website.herokuapp.com/MM/DNVS/ship_003.tar.gz"
	downloadSS "shipsec8" "https://suitesparse-collection-website.herokuapp.com/MM/DNVS/shipsec8.tar.gz"
	downloadSS "tmt_sym" "https://suitesparse-collection-website.herokuapp.com/MM/CEMW/tmt_sym.tar.gz"
	downloadSS "x104" "https://suitesparse-collection-website.herokuapp.com/MM/DNVS/x104.tar.gz"
	echo ""
	echo "Finished downloading."
	echo "The matrices are available in: ${ALLMATRICES_DIR}"
	echo ""
	if [[ ! -d "${SUITESPARSE_DIR}" ]]; then
		mkdir "${SUITESPARSE_DIR}" || exit 1
	fi
	for mat in $SUITESPARSE_SET; do
		if [ ! -f ${SUITESPARSE_DIR}/${mat}.mtx ]; then
			ln ${mat}.mtx ${SUITESPARSE_DIR}/${mat}.mtx || exit 1
		fi
	done
	echo "Finished SuiteSparse data set in ${SUITESPARSE_DIR}"
	if [[ ! -d "${ICHOL_DIR}" ]]; then
		mkdir "${ICHOL_DIR}" || exit 1
	fi
	for mat in $ICHOL_SET; do
		if [ ! -f ${mat}_postChol.mtx ]; then
			${OSP_DIR}/build/examples/post_incomplete_cholesky ${mat}.mtx || exit 1
		fi
		if [ ! -f ${ICHOL_DIR}/${mat}_postChol.mtx ]; then
			ln ${mat}_postChol.mtx ${ICHOL_DIR}/${mat}_postChol.mtx || exit 1
		fi
	done
	echo "Finished incomplete Cholesky data set in ${ICHOL_DIR}"
	if [[ ! -d "${ERDOS_RENYI_DIR}" ]]; then
		mkdir "${ERDOS_RENYI_DIR}" || exit 1
	fi
	cd "${ERDOS_RENYI_DIR}" || exit 1
	${OSP_DIR}/build/examples/gen_Erdos-Renyi_graph 100000 10 10 || exit 1
	${OSP_DIR}/build/examples/gen_Erdos-Renyi_graph 100000 50 10 || exit 1
	${OSP_DIR}/build/examples/gen_Erdos-Renyi_graph 100000 200 10 || exit 1
	echo "Finished Erdos--Renyi data set in ${ERDOS_RENYI_DIR}"
	if [[ ! -d "${NARROW_BANDWIDTH_DIR}" ]]; then
		mkdir "${NARROW_BANDWIDTH_DIR}" || exit 1
	fi
	cd "${NARROW_BANDWIDTH_DIR}" || exit 1
	${OSP_DIR}/build/examples/gen_near_diag_random_graph 100000 0.14 10 10 || exit 1
	${OSP_DIR}/build/examples/gen_near_diag_random_graph 100000 0.05 20 10 || exit 1
	${OSP_DIR}/build/examples/gen_near_diag_random_graph 100000 0.03 42 10 || exit 1
	echo "Finished narrow bandwidth data set in ${NARROW_BANDWIDTH_DIR}"
	if [[ ! -d "${METIS_DIR}" ]]; then
		mkdir "${METIS_DIR}" || exit 1
	fi
	if [[ ! "${SYMPILER_METIS_BIN}" = "" ]]; then
		cd "${ALLMATRICES_DIR}" || exit 1
		for mat in $METIS_SET; do
			if [ ! -f ${mat}_metis.mtx ]; then
				${SYMPILER_METIS_BIN} ${mat}.mtx || exit 1
			fi
			if [ ! -f ${METIS_DIR}/${mat}_metis.mtx ]; then
				ln ${mat}_metis.mtx ${METIS_DIR}/${mat}_metis.mtx || exit 1
			fi
		done
		echo "Finished Metis data set in ${METIS_DIR}"
	else
		echo "Add the absolute path to the Sympiler metis matrices generation binary to the script"
		exit 3
	fi
	exit 0
else
	echo "'yes' is required to continue."
	exit 2
fi


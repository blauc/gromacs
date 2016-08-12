/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include <map>
#include "emscatteringfactors.h"

namespace gmx
{

namespace externalpotential
{

real atomicnumber2emscatteringfactor (int atomic_number)
{
    static std::map < int, real > atomicnumber2factor {{
                                                           0, 1.0
                                                       }, {
                                                           1, 0.0444872
                                                       }, {
                                                           2, 0.0277892
                                                       }, {
                                                           3, 1.71605
                                                       }, {
                                                           4, 1.48231
                                                       }, {
                                                           5, 1.24151
                                                       }, {
                                                           6, 0.999998
                                                       }, {
                                                           7, 0.780615
                                                       }, {
                                                           8, 0.627736
                                                       }, {
                                                           9, 0.5166
                                                       }, {
                                                           10, 0.433072
                                                       }, {
                                                           11, 3.61627
                                                       }, {
                                                           12, 4.29754
                                                       }, {
                                                           13, 5.50428
                                                       }, {
                                                           14, 5.3601
                                                       }, {
                                                           15, 4.79207
                                                       }, {
                                                           16, 4.24146
                                                       }, {
                                                           17, 3.76566
                                                       }, {
                                                           18, 3.34392
                                                       }, {
                                                           19, 12.6911
                                                       }, {
                                                           20, 15.5702
                                                       }, {
                                                           21, 13.6654
                                                       }, {
                                                           22, 12.1117
                                                       }, {
                                                           23, 10.7941
                                                       }, {
                                                           24, 7.69918
                                                       }, {
                                                           25, 8.86039
                                                       }, {
                                                           26, 8.13567
                                                       }, {
                                                           27, 7.44624
                                                       }, {
                                                           28, 6.81748
                                                       }, {
                                                           29, 4.93804
                                                       }, {
                                                           30, 5.84716
                                                       }, {
                                                           31, 8.00058
                                                       }, {
                                                           32, 8.64189
                                                       }, {
                                                           33, 8.53914
                                                       }, {
                                                           34, 8.27235
                                                       }, {
                                                           35, 7.92845
                                                       }, {
                                                           36, 7.55695
                                                       }, {
                                                           37, 21.8987
                                                       }, {
                                                           38, 27.024
                                                       }, {
                                                           39, 25.4741
                                                       }, {
                                                           40, 23.5588
                                                       }, {
                                                           41, 18.2075
                                                       }, {
                                                           42, 16.7869
                                                       }, {
                                                           43, 18.6131
                                                       }, {
                                                           44, 14.4776
                                                       }, {
                                                           45, 13.5372
                                                       }, {
                                                           46, 9.12112
                                                       }, {
                                                           47, 11.9247
                                                       }, {
                                                           48, 13.5244
                                                       }, {
                                                           49, 17.235
                                                       }, {
                                                           50, 18.7
                                                       }, {
                                                           51, 19.157
                                                       }, {
                                                           52, 19.16
                                                       }, {
                                                           53, 18.8988
                                                       }, {
                                                           54, 18.4903
                                                       }, {
                                                           55, 42.5438
                                                       }, {
                                                           56, 52.3393
                                                       }, {
                                                           57, 50.1869
                                                       }, {
                                                           58, 47.8218
                                                       }, {
                                                           59, 45.5389
                                                       }, {
                                                           60, 43.1338
                                                       }, {
                                                           61, 41.5183
                                                       }, {
                                                           62, 39.7828
                                                       }, {
                                                           63, 38.1517
                                                       }, {
                                                           64, 37.0207
                                                       }, {
                                                           65, 35.2697
                                                       }, {
                                                           66, 33.898
                                                       }, {
                                                           67, 32.5976
                                                       }, {
                                                           68, 31.1725
                                                       }, {
                                                           69, 30.0192
                                                       }, {
                                                           70, 28.2632
                                                       }, {
                                                           71, 28.0326
                                                       }, {
                                                           72, 26.9922
                                                       }, {
                                                           73, 25.7361
                                                       }, {
                                                           74, 24.704
                                                       }, {
                                                           75, 23.8684
                                                       }, {
                                                           76, 22.6674
                                                       }, {
                                                           77, 21.4947
                                                       }, {
                                                           78, 18.5165
                                                       }, {
                                                           79, 17.6589
                                                       }, {
                                                           80, 18.9234
                                                       }, {
                                                           81, 25.6889
                                                       }, {
                                                           82, 26.9626
                                                       }, {
                                                           83, 26.9474
                                                       }, {
                                                           84, 23.8629
                                                       }, {
                                                           85, 28.734
                                                       }, {
                                                           86, 28.7999
                                                       }, {
                                                           87, 50.9563
                                                       }, {
                                                           88, 62.8004
                                                       }, {
                                                           89, 66.5727
                                                       }, {
                                                           90, 64.8004
                                                       }, {
                                                           91, 60.5933
                                                       }, {
                                                           92, 57.8371
                                                       }, {
                                                           93, 55.6476
                                                       }, {
                                                           94, 51.8924
                                                       }, {
                                                           95, 49.8216
                                                       }, {
                                                           96, 49.1348
                                                       }, {
                                                           97, 47.6216
                                                       }, {
                                                           98, 44.1356
                                                       }, {
                                                           99, 42.4524
                                                       }, {
                                                           100, 40.7625
                                                       }, {
                                                           101, 39.2219
                                                       }, {
                                                           102, 37.6153
                                                       }, {
                                                           103, 39.7949
                                                       }};
    return atomicnumber2factor[atomic_number];
}


} /* externalpotential */
} /* gmx */

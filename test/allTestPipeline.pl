#!/usr/bin/env perl

###################################################################################################################################
#
# Copyright 2014-2018 IRD-CIRAD-INRA-ADNid
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
# Version 1 written by Cecile Monat, Ayite Kougbeadjo, Christine Tranchant, Cedric Farcy, Mawusse Agbessi, Maryline Summo, and Francois Sabot
# Version 2 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Enrique Ortega-Abboud, Julie Orjuela-Bouniol, Sebastien Ravel, Souhila Amanzougarene, and Francois Sabot
# Version 3 written by Cecile Monat, Christine Tranchant, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
#
###################################################################################################################################

use strict;
use warnings;
use localConfig;


#####################
## PIPELINE TEST
#####################

print "\n\n#################################################\n";
print "#### PIPELINE TEST \n";
print "#################################################\n";

system("perl $toggle/test/pipelines/pairedOneIndividuArcad.pl") and warn "ERROR: $0: Cannot run test for pairedOneIndividuArcad.pl  \n$!\n";
system("perl $toggle/test/pipelines/pairedTwoIndividusGzippedIrigin.pl") and warn "ERROR: $0: Cannot run test for pairedTwoIndividusGzippedIrigin.pl  \n$!\n";
system("perl $toggle/test/pipelines/pairedTwoIndividusIrigin.pl") and warn "ERROR: $0: Cannot run test for pairedTwoIndividusIrigin.pl  \n$!\n";
system("perl $toggle/test/pipelines/pairedTwoIndividusIriginSGE.pl") and warn "ERROR: $0: Cannot run test for pairedTwoIndividusIriginSGE.pl  \n$!\n";
system("perl $toggle/test/pipelines/processRadtags.pl") and warn "ERROR: $0: Cannot run test for processRadtags.pl  \n$!\n";
system("perl $toggle/test/pipelines/rnaSeqPairedOneIndividu.pl") and warn "ERROR: $0: Cannot run test for rnaSeqPairedOneIndividu.pl  \n$!\n";
system("perl $toggle/test/pipelines/samBam.pl") and warn "ERROR: $0: Cannot run test for samBam.pl  \n$!\n";
system("perl $toggle/test/pipelines/singleOneIndividuIrigin.pl") and warn "ERROR: $0: Cannot run test for singleOneIndividuIrigin.pl  \n$!\n";
system("perl $toggle/test/pipelines/singleTwoIndividuIrigin.pl") and warn "ERROR: $0: Cannot run test for singleTwoIndividuIrigin.pl  \n$!\n";


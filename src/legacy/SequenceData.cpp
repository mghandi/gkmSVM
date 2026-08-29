/* SequenceData.cpp
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "SequenceData.h"
#include "global.h"


CSequenceData::CSequenceData(void)
{

	
	seq=NULL; 
	iseq=NULL; 

	nchrs=0; // number of chromosomes 
	maxchrsize=0; 
	chrsize=NULL; 

}


CSequenceData::~CSequenceData(void)
{
}

void CSequenceData::readFrombpmap(char *bpmapfn)
{
}
void CSequenceData::readFromFASTA(char *genomefn)
{
}

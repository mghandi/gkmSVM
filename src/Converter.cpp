/* Converter.cpp: implementation of the CConverter class.
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

//#include "stdafx.h"
#include "Converter.h"
#include "global.h"
#include <string>
#include <string.h>
#include <ctype.h>
//#include <stdio.h>
//#include <ctype.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CConverter::CConverter()
{
  b=4;
  explicitComplement=0;
  alphabet[0]='A';
  alphabet[1]='C';
  alphabet[2]='G';
  alphabet[3]='T';
  
  init();
  
}

void CConverter::init(){
  
  //char ci;
  int ici;
  //for (ci='A';ci<'z';ci++) cidx[ci]=0;
  for (ici=0;ici<256;ici++) cidx[ici]=0;
  
  //	cidx['a'] = 0; 		cidx['c'] = 1; 		cidx['g'] = 2; 		cidx['t'] = 3;
  //	cidx['A'] = 0; 		cidx['C'] = 1; 		cidx['G'] = 2; 		cidx['T'] = 3;
  for(ici=0;ici<b;ici++){
    cidx[toupper(alphabet[ici])]=ici;
    cidx[tolower(alphabet[ici])]=ici;
  }
  
  icidx = new char[b];
  icidxL = new char[b];
  
  //	icidx[0]  = 'A';		icidx[1]  = 'C';		icidx[2]  = 'G';		icidx[3]  = 'T';
  //	icidxL[0] = 'a';		icidxL[1] = 'c';		icidxL[2] = 'g';		icidxL[3] = 't';
  for(ici=0;ici<b;ici++){
    icidx[ici]=toupper(alphabet[ici]);
    icidxL[ici]=tolower(alphabet[ici]);
    
    if (!explicitComplement) {
      bidcompl[ici]=b-ici-1; // the built-in rule: DNA/RNA in ACGT/ACGU order; any other 4-letter alphabet pairs 1st<->4th, 2nd<->3rd
      if (b==16){
        // special case for dinucleotides
        bidcompl[ici]= ((3-(ici&3))<<2) + (3-((ici&12)>>2)) ; //only good for DNA
      }
    }
  }
  hasComplement = explicitComplement || (b==4) || (b==16);
  
  //	bidcompl[0] = 3;
  //	bidcompl[1] = 2;
  //	bidcompl[2] = 1;
  //	bidcompl[3] = 0;
  
  for (ici=0;ici<256;ici++) //(ci='A';ci<'z';ci++)
  {
    //bcompl[ci] = icidx[bidcompl[cidx[ci]]];
    bcompl[ici] = icidx[bidcompl[cidx[ici]]];
  }
  
  for (ici=0;ici<256;ici++)
  {
    isACGT[ici] = 0;
    isInAlphabet[ici] = 0;
  }
  isACGT['a'] = 1; 		isACGT['c'] = 1; 		isACGT['g'] = 1; 		isACGT['t'] = 1;
  isACGT['A'] = 1; 		isACGT['C'] = 1; 		isACGT['G'] = 1; 		isACGT['T'] = 1;
  
  for(ici=0;ici<b;ici++){
    isInAlphabet[toupper(alphabet[ici])]=1;
    isInAlphabet[tolower(alphabet[ici])]=1;
  }
  
  
}
CConverter::~CConverter()
{
  delete []icidx; 
  delete []icidxL; 
  
}


/*

aa	0
ca	1
ga	2
ta	3
ac	4
cc	5
gc	6
tc	7
ag	8
cg	9
gg	10
tg	11
at	12
ct	13
gt	14
tt	15
*/


dinuclId CConverter::dnidx(const char *dn) const
{
  return cidx[(unsigned char)dn[0]]+b*cidx[(unsigned char)dn[1]];
}




void CConverter::convertBasetoDinucl(baseId x[], dinuclId y[], int N) // x[0..N], y[0..N-1]
{
  for(int i=0;i<N;i++)
    y[i] = x[i]+b*x[i+1];
}

void CConverter::convertBasetoDinucl(char x[], dinuclId y[], int N) // x[0..N], y[0..N-1]
{
  for(int i=0;i<N;i++)
    y[i] = cidx[(unsigned char)x[i]]+b*cidx[(unsigned char)x[i+1]];
}


void CConverter::resetToDNA(){
  b=4;
  explicitComplement=0;
  alphabet[0]='A';
  alphabet[1]='C';
  alphabet[2]='G';
  alphabet[3]='T';
  delete []icidx;
  delete []icidxL;
  init();
}

int CConverter::readAlphabetFile(const char *FN, int MAX_ALPHABET_SIZE_copy){
  int named = setNamedAlphabet(FN);
  if (named >= 0) { snprintf(globtmpstr, GKM_TMPSTR_LEN,"Alphabet: %s (%d symbols)\n", FN, b);Printf(globtmpstr); return named; }
  FILE *f= fopen(FN,"r");
  if (f==NULL){
    snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n ERROR: cannot open alphabet file %s\n", FN);Printf(globtmpstr);
    return 1;
  }
  // One symbol per line; an optional second symbol on the line is its complement ("A T"). When
  // every symbol has a complement, reverse complements can be added (addRC); otherwise addRC is
  // turned off by the caller. Lines starting with '#' are comments.
  char sline[1000+3];
  int nb=0;
  char partner[256]; int npairs=0;
  while(fgets(sline, 1000, f)!=NULL){
    if (sline[0]=='\n' || sline[0]=='\r' || sline[0]==0 || sline[0]=='#' || sline[0]==' ' || sline[0]=='\t') continue;
    if (nb>=MAX_ALPHABET_SIZE_copy) { nb++; break; }
    alphabet[nb]=sline[0];
    partner[nb]=0;
    int i=1; while (sline[i]==' '||sline[i]=='\t') i++;
    if (sline[i]!=0 && sline[i]!='\n' && sline[i]!='\r') { partner[nb]=sline[i]; npairs++; }
    nb++;
  }
  fclose(f);
  snprintf(globtmpstr, GKM_TMPSTR_LEN,"Alphabet Size = %d\n",nb);Printf(globtmpstr);
  if(nb>MAX_ALPHABET_SIZE_copy){
    snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n ERROR: alphabet size (>%d) is greater than the supported maximum of %d symbols.\n", MAX_ALPHABET_SIZE_copy, MAX_ALPHABET_SIZE_copy);Printf(globtmpstr);
    return 1;
  }
  if(nb<2){
    Printf("\n ERROR: alphabet file must contain at least two symbols, one per line.\n");
    return 1;
  }
  for (int i=0;i<nb;i++) for (int j=0;j<i;j++) if (toupper(alphabet[i])==toupper(alphabet[j])) {
    snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n ERROR: symbol '%c' appears twice in the alphabet file.\n", alphabet[i]);Printf(globtmpstr);
    return 1;
  }
  return setAlphabet(alphabet, nb, npairs ? partner : NULL);
}

int CConverter::setAlphabet(const char *symbols, int nb, const char *partner){
  if (symbols != alphabet) memcpy(alphabet, symbols, nb);
  explicitComplement = 0;
  if (partner != NULL) {
    // every symbol must have a partner that is itself in the alphabet, and pairing must be symmetric
    for (int i=0;i<nb;i++) {
      if (partner[i]==0) { snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n ERROR: symbol '%c' has no complement while others do; declare a complement for every symbol or for none.\n", alphabet[i]);Printf(globtmpstr); return 1; }
      int j; for (j=0;j<nb;j++) if (toupper(alphabet[j])==toupper(partner[i])) break;
      if (j==nb) { snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n ERROR: complement '%c' of '%c' is not in the alphabet.\n", partner[i], alphabet[i]);Printf(globtmpstr); return 1; }
      bidcompl[i]=j;
    }
    for (int i=0;i<nb;i++) if (bidcompl[bidcompl[i]]!=i) { snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n ERROR: complement pairs are not symmetric ('%c' -> '%c' -> '%c').\n", alphabet[i], alphabet[bidcompl[i]], alphabet[bidcompl[bidcompl[i]]]);Printf(globtmpstr); return 1; }
    explicitComplement = 1;
  }
  b=nb;
  delete []icidx;
  delete []icidxL;
  init();
  return 0;
}

// Built-in alphabets, selectable by keyword instead of a file: "dna", "rna", "protein".
int CConverter::setNamedAlphabet(const char *name){
  std::string n(name); for (size_t i=0;i<n.size();i++) n[i]=tolower(n[i]);
  if (n=="dna") { resetToDNA(); return 0; }
  if (n=="rna") { char rna[4]={'A','C','G','U'}; char part[4]={'U','G','C','A'}; return setAlphabet(rna, 4, part); }
  if (n=="protein") { const char *aa="ACDEFGHIKLMNPQRSTVWY"; return setAlphabet(aa, 20, NULL); }
  return -1; // not a keyword
}


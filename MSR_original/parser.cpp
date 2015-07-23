
/*
 *  parser.cpp
 *
 *
 *  Created by Antoni Bueades Capó.
 *
 *  Copyright 2012 IPOL Image Processing On Line http://www.ipol.im/
 *
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

#include "parser.h"

void  printusage(char *pname,
                 char *gp,
                 std::vector<OptStruct*>  &opt,
                 std::vector<ParStruct*>  &par)
{
    
    
    int npar = par.size();
    
    //! print function prototype
    //! "usage:" function_name [-r r] [-g g] [-b b]  par1 par2 ... parn
    printf("\nusage: %s ", pname);
    for(int i=0; i < (int) strlen(gp); i++)
        if (gp[i] != ':')
        {
            printf("[-%c",gp[i]);
            
            if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':') printf(" %c] ", gp[i]);
            else printf("] ");
            
        }
    
    for(int i=0; i < npar; i++)
        printf(" %s ", par[i]->name);
    
    printf("\n");
    
    
    //! print options with associated descriptions and defaulted values
    int j=0;
    for(int i=0; i < (int) strlen(gp); i++)
        if (gp[i] != ':')
        {
            printf("\t-%c",gp[i]);
            
            if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':') {
                
                printf("  %c\t %s ", gp[i], opt[j]->comment);
                if (opt[j]->defvalue != NULL) printf("(Default: %s)",opt[j]->defvalue);
                
                printf("\n");
                
            }
            else printf("\t %s \n", opt[j]->comment);
            
            j++;
        }
    
    //! print mandatory parameters with associated descriptions
    for(int i=0; i < npar; i++)
    {
        printf("\t%s",par[i]->name);
        printf("\t %s\n", par[i]->comment);
    }
    
    
    
}




int parsecmdline(char *pname,
                 char *function,
                 int argc, char **argv,
                 std::vector <OptStruct*> & opt,
                 std::vector <ParStruct*> & par)
{
    
    //! number of options and obligatory parameters
    int nopt = opt.size();
    int npar = par.size();
    
    
    //! size of auxiliar parameter needed by getopt
    //! the length is at maximum 2*nopt since each option contains one caracter
    //! or two if the option requires a values
    char *gp = new char[2*nopt+1];
    gp[0]='\0';
    
    
    //! clear options and concatenate option identifiers into gp vector
    for(int i=0; i < nopt; i++) { opt[i]->flag = 0; opt[i]->value=NULL; strcat(gp, opt[i]->gp);}
    
    //! clear necessary parameter values
    for(int i=0; i < npar; i++) { par[i]->value = NULL;}
    
    
    //! in this way getopt doesn't print any information and errors are
    //! treated and printed by this program. 
    opterr = 0;
    
    
    //! getopt return the following option of the concatenated vector gp
    //! or -1 if all options have already been treated
    //! it also fills "optarg" with the associated value to this option as it is given in argv 
    int c;
    while ((c = getopt (argc, argv, gp)) != -1)
    {
        
        //! set current option given by getopt
        //! if getopt finds an option which was not in the option of the program (vector gp)
        //! returns ?
        int j=0;
        for(unsigned int i=0; i < strlen(gp); i++)
            if (c == gp[i])
            {
                opt[j]->flag = 1;
                opt[j]->value = optarg;  // getopt fills "optarg" with the associated value to this option as it is given in argv
                
                break;
                
            } else if (gp[i] != ':') j++;
        
        
        //! option found in console command was not one of the parameters of our program
        //! or should have a mandatory values which is not provided
        if (c == '?')
        {
            
            
            //! when getopt encounters an unknown option character or an option with a missing required argument
            //! it stores that option character in optopt variable
            unsigned int i = 0;
            for(i=0; i < strlen(gp); i++)
                if (optopt == gp[i])
                {
                    printf("\n%s: %s\n", pname, function);
                    fprintf (stderr, "\nerror: option -%c requires an argument.\n", optopt);
                    break;
                }
            
            if (i == strlen(gp)) {
                printf("\n%s: %s\n", pname, function);
                fprintf (stderr, "\nerror: unknown option `-%c'.\n", optopt);
            }
            
            printusage(pname, gp,  opt,  par);
            delete[] gp;
            return 0;
            
        }
        
        
        
    }
    
    
    //! Setting default values for non selected options
    for(int j=0; j < nopt; j++)
        if (opt[j]->flag == 0 && opt[j]->defvalue != NULL) opt[j]->value =  opt[j]->defvalue;
    
    
    //! Check remaining words in command after reading option
    if (argc - optind != npar)
    {
        printf("\n%s: %s\n", pname, function);
        fprintf (stderr, "\nerror: incorrect number of parameters\n");
        printusage(pname, gp,  opt,par);
        delete[] gp;
        return 0;
    }
    
    //! Read mandatory parameter values
    int i=0;
    for (int index = optind; index < argc ; index++, i++){
        par[i]->value = argv[index];
    }
    
    
    delete[] gp;
    return 1;
    
    
}


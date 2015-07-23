/*
 *  parser.h
 *  
 *
 *  Created by Antoni Buades 2012.
 *
 * Copyright 2012 IPOL Image Processing On Line http://www.ipol.im/
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

#ifndef _PARSER_H
#define _PARSER_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <vector>
#include <getopt.h>


//! structure for parameters and options which are
//! optional in the code or they already have a default value
typedef struct optstruct
{
    char *gp;           //! string of two letters  "a:"  as necessary for using getopt afterwards
    //! the ":" indicates that the activation of the option requires a value for the parameter
    //! and "a" that this option is activated by "-a " in the command
    
    int flag;           //! flag indicating that the option has been activated
    
    char *defvalue;     //! default value for this parameter if not modified by console
    
    char *value;        //! value of the associated parameter to current option
    
    char *comment;      //! comment that appears by console
    
} OptStruct;



//! structure for necessary parameters of the method
typedef struct parstruct
{
    char * name;
    char * value;       //! value of the parameter
    char * comment;     //! comment that appears by console
    
} ParStruct;



int parsecmdline(char *pname,
                 char *function,
                 int argc, char **argv,
                 std::vector <OptStruct*> & opt,
                 std::vector <ParStruct*> & par);



#endif /* _PARSER_H */





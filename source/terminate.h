/******************************************************************************/
/**                                                                          **/
/**     Program Termination                                                  **/
/**             functions to print out warning messages,                     **/
/**             and to terminate code execution.                             **/
/**                                                                          **/
/******************************************************************************/

//#include <ostream>
#include <sstream>

#ifndef CENS_TOPLEVEL
extern ostringstream message;
#endif

/**************************************/
/*      cens.cpp                      */
/**************************************/
void    WarningMessage     ();
void    Notice             (std::string);
int     TerminateCode      (std::string);

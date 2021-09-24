#pragma once

#include <iostream>
#include <string>
#include <experimental/filesystem>



using namespace std;

template< typename T>
struct base_options
{
    string outfile = "ecgsyn.dataecg";/*  Output data file                   */
    int N = 256;                   /*  Number of heart beats              */
    int sfecg = 256;               /*  ECG sampling frequency             */
    int sf = 256;                  /*  Internal sampling frequency        */
    double Anoise = 0.0;           /*  Amplitude of additive uniform noise*/
    double hrmean = 60.0;          /*  Heart rate mean                    */
    double hrstd = 1.0;            /*  Heart rate std                     */
    double flo = 0.1;              /*  Low frequency                      */
    double fhi = 0.25;             /*  High frequency                     */
    double flostd = 0.01;          /*  Low frequency std                  */
    double fhistd = 0.01;          /*  High frequency std                 */
    double lfhfratio = 0.5;        /*  LF/HF ratio                        */

    int Necg = 0;                  /*  Number of ECG outputs              */
    int mstate = 3;                /*  System state space dimension       */
    double xinitial = 1.0;         /*  Initial x co-ordinate value        */
    double yinitial = 0.0;         /*  Initial y co-ordinate value        */
    double zinitial = 0.04;        /*  Initial z co-ordinate value        */
    int seed = 1;                  /*  Seed                               */
    long rseed;

    base_options ()
    {

    }


};



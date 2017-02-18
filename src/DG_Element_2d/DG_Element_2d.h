/**
 * @file DG_Element_2d.h
 * @Synopsis  This is the header file for the class which defines the DG_Element_2d class. This class contains all the
 * declarations which are necessary for defining a general 2d element. This class defines all the necessities for a
 * structured 2d DG element.
 * @author Kaushik Kulkarni
 * @version 1.0
 * @date 2017-02-18
 */



/* Copyright (c) 2016, Shivasubramanian Gopalakrishnan and Kaushik Kulkarni.
  All rights reserved.

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

     -Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
     -Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS  **AS IS** AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     */

#ifndef DG_ELEMENT_2D_H
#define DG_ELEMENT_2D_H

#include <map>
#include <string>

using namespace std;


class DG_Element_2d {
private:
    
    int N; /// This represents the order of polynomial interpolation which is done while solving the DG equations.
    float x_start, x_end, y_start, y_end; /// Since, this is a structured rectangular element. These variables define the position and size of the element.

public:

    float *X, *Y;
    map<string, float*> variable; /// This is the map which would contain all the variable that are needed by the problem.


    map<string, float**> boundaryTop;
    map<string, float**> boundaryRight;
    map<string, float**> boundaryLeft;
    map<string, float**> boundaryBottom;
    
    map<string, float**> neighboringTop;
    map<string, float**> neighboringRight;
    map<string, float**> neighboringBottom;
    map<string, float**> neighboringLeft;

    DG_Element_2d(int _N, float x1, float y1, float x2, float y2);
    void addVariable_withBoundary(string v);
    void addVariable_withoutBoundary(string v);

};

#endif

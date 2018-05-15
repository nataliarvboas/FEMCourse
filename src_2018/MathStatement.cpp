/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "MathStatement.h"

MathStatement::MathStatement() {
}

MathStatement::MathStatement(const MathStatement &copy) {
    matid = copy.matid;
    nstate = copy.nstate;
}

MathStatement &MathStatement::operator=(const MathStatement &copy) {
    matid = copy.matid;
    nstate = copy.nstate;
    return *this;
}

MathStatement::~MathStatement() {
}


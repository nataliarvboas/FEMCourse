/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <stdio.h>
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "GeoElementSide.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "Topology1d.h"
#include "TopologyQuad.h"
#include "TopologyTriangle.h"
#include "TopologyTetrahedron.h"
#include <algorithm>


GeoElementSide::GeoElementSide() {

}

GeoElementSide::GeoElementSide(const GeoElementSide &copy) {
    fElement = copy.fElement;
    fSide = copy.fSide;
}

GeoElementSide &GeoElementSide::operator=(const GeoElementSide &copy) {
    fElement = copy.fElement;
    fSide = copy.fSide;
    return *this;
}

GeoElementSide GeoElementSide::Neighbour() const {
    return fElement ? fElement->Neighbour(fSide) : GeoElementSide();
}

void GeoElementSide::SetNeighbour(const GeoElementSide &neighbour) {
    fElement->SetNeighbour(fSide, neighbour);
}

bool GeoElementSide::IsNeighbour(const GeoElementSide &candidate) {
    if (candidate == *this) return true;
    GeoElementSide neighbour = Neighbour();
    if (!(neighbour.Element() != 0 && neighbour.Side() > -1)) return false;
    while ((neighbour == *this) == false) {
        if (candidate == neighbour) return true;
        neighbour = neighbour.Neighbour();
    }
    return false;
}

void GeoElementSide::InsertConnectivity(GeoElementSide &candidate) {
    if (this->IsNeighbour(candidate)) return;
    GeoElementSide neigh1 = Neighbour();
    GeoElementSide neigh2 = candidate.Neighbour();
    neigh1.SetNeighbour(neigh2);
    candidate.SetNeighbour(neigh1);
}

void GeoElementSide::AllNeighbours(std::vector<GeoElementSide> &allneigh) {
    GeoElementSide neigh = Neighbour();

    while (neigh != *this) {
        allneigh.push_back(neigh);
        neigh = neigh.Neighbour();
    }
}

void GeoElementSide::Intersect(const std::vector<int> &one, const std::vector<int> &two, std::vector<int> &result) {
    int firstc, secondc, nfirst, nsecond;
    nfirst = one.size();
    nsecond = two.size();
    firstc = 0;
    secondc = 0;
    while (firstc < nfirst && secondc < nsecond) {
        while (firstc < nfirst && one[firstc] < two[secondc]) {
            firstc++;
        }
        if (firstc == nfirst) break;
        while (secondc < nsecond && two[secondc] < one[firstc]) {
            secondc++;
        }
        if (firstc < nfirst && secondc < nsecond && one[firstc] == two[secondc]) {
            result.push_back(one[firstc]);
            firstc++;
            secondc++;
        }
    }
}

void GeoElementSide::Intersect(const std::vector<int> &one, const std::vector<int> &two, const std::vector<int> &three, std::vector<int> &result) {
    int firstc, secondc, thirdc, nfirst, nsecond, nthird;
    nfirst = one.size();
    nsecond = two.size();
    nthird = three.size();
    firstc = 0;
    secondc = 0;
    thirdc = 0;
    while (firstc < nfirst && secondc < nsecond && thirdc < nthird) {
        while (firstc < nfirst && (one[firstc] < two[secondc] || one[firstc] < three[thirdc])) {
            firstc++;
        }
        if (firstc == nfirst)break;
        while (secondc < nsecond && (two[secondc] < one[firstc] || two[secondc] < three[thirdc])) {
            secondc++;
        }
        if (secondc == nsecond) break;
        while (thirdc < nthird && (three[thirdc] < one[firstc] || three[thirdc] < two[secondc])) {
            thirdc++;
        }
        if (firstc < nfirst && secondc < nsecond && thirdc < nthird && one[firstc] == two[secondc] && one[firstc] == three[thirdc]) {
            result.push_back(one[firstc]);
            firstc++;
            secondc++;
            thirdc++;
        }
    }
}

void GeoElementSide::ComputeNeighbours(std::vector<GeoElementSide> &compneigh) {
    if (fSide < fElement->NCornerNodes()) {
        AllNeighbours(compneigh);
        return;
    }
    int nsnodes = fElement->NSideNodes(fSide);
    std::vector<GeoElementSide> GeoElSideSet;
    std::vector<std::vector<int>> GeoElSet;
    GeoElSet.resize(27);
    int in;
    VecInt nodeindexes(nsnodes);
    for (in = 0; in < nsnodes; in++) {
        nodeindexes[in] = fElement->NodeIndex(fElement->SideNodeIndex(fSide, in));       
        int locnod = fElement->SideNodeIndex(fSide, in);
        GeoElSideSet.resize(0);
        GeoElementSide locside(fElement, locnod);
        locside.AllNeighbours(GeoElSideSet);
        int nel = GeoElSideSet.size();
        int el;
        for (el = 0; el < nel; el++) {
            GeoElSet[in].push_back(GeoElSideSet[el].Element()->GetIndex());
        }
        std::sort(GeoElSet[in].begin(), GeoElSet[in].end());
    }
    std::vector<int> result;
    switch (nsnodes) {
        case 1:
        {
            result = GeoElSet[0];
        }
            break;
        case 2:
            Intersect(GeoElSet[0], GeoElSet[1], result);
            break;
        case 3:
            Intersect(GeoElSet[0], GeoElSet[1], GeoElSet[2], result);
            break;
        case 4:
        {
            std::vector<int> inter1, inter2;
            Intersect(GeoElSet[0], GeoElSet[2], inter1);
            if (inter1.size() == 0) break;
            Intersect(GeoElSet[1], GeoElSet[3], inter2);
            if (inter2.size() == 0) break;
            Intersect(inter1, inter2, result);
        }
            break;
        default:
        {
            std::vector<int> inter1, inter2;
            inter1 = GeoElSet[0];
            for (in = 0; in < nsnodes - 1; in++) {
                inter2.resize(0);
                Intersect(inter1, GeoElSet[in + 1], inter2);
                if (inter2.size() == 0) break;
                inter1 = inter2;
            }
            result = inter2;
        }
    }
    int el, nel = result.size();
    GeoMesh * geoMesh = fElement->GetMesh();
    for (el = 0; el < nel; el++) {
        GeoElement * gelResult = geoMesh->Element(result[el]);
        int whichSd = gelResult->WhichSide(nodeindexes);
        if (whichSd > 0) {
            compneigh.push_back(GeoElementSide(gelResult, whichSd));
        }
    }

}

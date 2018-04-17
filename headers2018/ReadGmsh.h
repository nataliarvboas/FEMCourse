//
//  ReadGmsh.h
//  FemSC
//
//  Created by Philippe Devloo on 17/04/18.
//

#ifndef ReadGmsh_h
#define ReadGmsh_h

class ReadGmsh
{
    
    /// reads the mesh contained in the file and fill the geometric mesh
    static void Read(GeoMesh &gmesh, std::string &filename);
    
};
#endif /* ReadGmsh_h */

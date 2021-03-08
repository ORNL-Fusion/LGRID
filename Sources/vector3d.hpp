//
//  vector3d.hpp
//  lgrid
//
//  Created by Cianciosa, Mark R. on 2/7/16.
//  Copyright © 2016 Cianciosa, Mark R. All rights reserved.
//

#ifndef vector3d_hpp
#define vector3d_hpp

#include <stdio.h>

//------------------------------------------------------------------------------
///  @brief A vector.
///
///  Class representing a vector in cartesian coordinates.
//------------------------------------------------------------------------------
class vector3d {
public:
/// x component.
    const double x;
/// y component.
    const double y;
/// z component.
    const double z;
    
    vector3d(const double x_, const double y_, const double z_);
    
    const double length() const;
    const vector3d normalize() const;
    
    static const vector3d cross_product(const vector3d &v1, const vector3d &v2);
    static const double dot_product(const vector3d &v1, const vector3d &v2);
};

const vector3d operator+(const vector3d &v1, const vector3d &v2);
const vector3d operator-(const vector3d &v1, const vector3d &v2);
const vector3d operator/(const vector3d &v1, const double d);

#endif /* vector3d_hpp */

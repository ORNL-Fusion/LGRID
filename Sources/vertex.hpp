//
//  vertex.h
//  lgrid
//
//  Created by Cianciosa, Mark R. on 10/12/15.
//  Copyright (c) 2015 Cianciosa, Mark R. All rights reserved.
//

#ifndef vertex_hpp
#define vertex_hpp

#include "vector3d.hpp"
#include <vector>

//------------------------------------------------------------------------------
///  @brief A vertex.
///
///  Class representing a vertex. Verticies are doublely linked lists containing
///  the current vector position.
//------------------------------------------------------------------------------
class vertex {
/// Vertex point.
    const vector3d point;

/// Vertex next vertex.
    vertex *next;
/// Vertex last vertex.
    vertex *last;
    
    public:
    
    vertex(const vector3d point_);
    ~vertex();
    
    void insert(vertex *v);
    
    static const vector3d get_normal(const vertex &v1, const vertex &v2);
    static const vector3d get_normal(const vertex &v1, const vertex &v2, const vertex &v3);
    
    static const double a(const vertex &v1, const vertex &v2);
    static const double b(const vertex &v1, const vertex &v2);
    static const double c(const vertex &v1, const vertex &v2);
    
    static const double x_min(const vertex &v1, const vertex &v2, const vector3d &p);
    static const double y_min(const vertex &v1, const vertex &v2, const vector3d &p);
    
    static const double length(const vertex &v1, const vertex &v2, const vector3d &p);
    static const double length(const vertex &v, const vector3d &p);
    
    static const bool in_range(const vertex &v1, const vertex &v2, const double x, const double y);

    static const double direction(const vertex &v, const vector3d &n, const vector3d &p);
    
    static void distance(vertex *start, const vector3d &p, std::vector<double> &d);
    
    static void print(vertex *start);
};

#endif /* defined(vertex_hpp) */

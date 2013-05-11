//
//  boid.h
//  
//
//  Created by Mark Whelan on 5/11/13.
//
//

#ifndef ____boid__
#define ____boid__

#include <iostream>

void GenerateBoids(R3Scene *scene, int quantity, double distAway);
void UpdateBoids(R3Scene *scene, double delta_time);

#endif /* defined(____boid__) */

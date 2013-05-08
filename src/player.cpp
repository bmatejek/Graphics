//
//  player.cpp
//  
//
//  Created by Mark Whelan on 5/7/13.
//
//

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "raytrace.h"
#include "player.h"
using namespace std;
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#endif
#define PI 3.141592655359
#define GRAV_CONSTANT 6.67428e-11
#define ADAPTIVE_THRESHOLD 1e-2
#define eps 2e-12



void UpdatePlayers(R3Scene *scene, double current_time, double delta_time, int integration_type) {
    
    scene->players[0]->pos += delta_time * (scene->players[0]->velocity * scene->players[0]->nose);
    
}
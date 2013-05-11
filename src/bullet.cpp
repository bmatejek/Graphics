//
//  bullet.cpp
//  
//
//  Created by Ethan Leeman on 5/7/13.
//
//

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "raytrace.h"
#include "bullet.h"
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



void UpdateBullets(R3Scene *scene, double current_time, double delta_time, int integration_type) {
    
    std::vector<int> tobedeleted;
    int i;
    for (i = 0; i < (int)scene->bullets.size(); i++) {
        R3Bullet *bullet = scene->bullets[i];
    
        bullet->position += delta_time * bullet->velocity;
        
        if (bullet->lifetimeactive) {
            
            bullet->lifetime -= delta_time;
            
            if (bullet->lifetime < 0) {
                
                tobedeleted.push_back(i);
                
            }
            
        }
        
    }
    for (unsigned int i = 0; i < tobedeleted.size(); i++) {
        //swap contents of partcile vector with last element
        R3Bullet *temp = scene->bullets.back();
        scene->bullets[scene->bullets.size() - 1] = scene->bullets[tobedeleted[i]];
        scene->bullets[tobedeleted[i]] = temp;
        
        // delete last element
        scene->bullets.pop_back();
    }
}

void RenderBullets(R3Scene *scene, double current_time, double delta_time)
{
    // Draw every particle
    
    // REPLACE CODE HERE
    glDisable(GL_LIGHTING);
    glPointSize(5);
    glBegin(GL_POINTS);

    
    for (int i = 0; i < (int)scene->bullets.size(); i += 10) {
        R3Bullet *bullet = scene->bullets[i];
        glColor3d(bullet->material->kd[0], bullet->material->kd[1], bullet->material->kd[2]);
        const R3Point& position = bullet->position;
        glVertex3d(position[0], position[1], position[2]);
    }
    glEnd();
}
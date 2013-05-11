//
//  bullet.h
//  
//
//  Created by Ethan Leeman on 5/7/13.
//
//

#ifndef _bullet_h
#define _bullet_h

void ShootBullet(R3Scene *scene);
void UpdateBullets(R3Scene *scene, double current_time, double delta_time, int integration_type);
void RenderBullets(R3Scene *scene, double current_time, double delta_time);

#endif

// Include file for ray tracing code

struct R3Intersect {
    bool intersected;
    R3Node *node;
    R3Point pos;
    R3Vector norm;
    double t;
};

R3Ray ConstructRayThroughPixel(int x, int y, int width, int height, R3Point eye, R3Vector towards, R3Vector up, R3Vector right, double xfov, double yfov);


R3Intersect PlaneIntersect(R3Ray r, R3Plane p);
R3Intersect BoxIntersect(R3Ray r, R3Box b) ;

// assumes mesh face is triangle
R3Intersect TriIntersect(R3Ray r, R3MeshFace f); 
R3Intersect TriMeshIntersect(R3Ray r, R3Mesh m);

double LightDist(R3Light *light, R3Point p);
R3Intersect SphereInt(R3Ray r, R3Sphere s);
R3Rgb Intensity(R3Light light, R3Point p);
R3Rgb DiffuseIntensity(R3Intersect sect, R3Light light, R3Material mat);
R3Rgb SpecularIntensity(R3Intersect sect, R3Light light, R3Material mat, R3Point camera);
R3Rgb ComputeRadiance(R3Scene *scene, R3Ray r, R3Intersect hit, int max_depth, double index);
R3Rgb ComputeRadiance(R3Scene *scene, R3Ray ray, int max_depth, double index);
R3Rgb Phong(R3Scene *scene, R3Ray r, R3Intersect hit);
R3Ray SpecularRay(R3Ray r, R3Intersect hit);
R3Ray TransmissionRay(R3Ray r, R3Intersect hit);
R3Ray RefractionRay(R3Ray r, R3Intersect hit, double indexinc, double indexout) ;
R3Vector LightVect(R3Light *light, R3Point p);
R3Ray LightRay(R3Light *light, R3Point p);


R3Intersect ComputeIntersect(R3Scene *scene, R3Node *node, R3Ray *r);
R3Rgb GetColor(R3Scene *scene, R3Ray r);
R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth,
                     int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection);
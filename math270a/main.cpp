#include <GL/glut.h> 
#include "ALGEBRA.h"
#include "GEOMETRY.h"
#include "DEFORMABLE_OBJECTS.h"
#include "FIXED_COROTATED.h"
#include <fstream>
#include <iostream>

using namespace ALGEBRA;
using namespace GEOMETRY;
using namespace std;

#define PI 3.14159265

typedef double T;
typedef VECTOR_3D<T> TV;
int winWidth,winHeight;
T costheta,sintheta;
MATRIX_3X3<T> ACCUMULATED=MATRIX_3X3<T>::Identity();
TV tsl(0.0, 0.0, 0.0);   //original eye position
const TV def_eye(0.0, 0.0, 5.0);   //original eye position
const TV def_ori(0.0, 1.0, 0.0);   //original orientation of camera
TV eye(def_eye);   //original eye position
TV ori(def_ori);   //original orientation of camera
TV u(0.0, 0.0, 5.0);     //current position
TV v1(0.0, 0.0, 5.0);    //last position
T one_over_sqrt_3=(T)1/sqrt((T)3);
T rotation_scale=(T)20;
T zoom_scale=(T)0.2;
T shift_scale=(T)0.05;

GLfloat light_diffuse[] = {0.0, 1.0, 0.0, 1.0};  /* Diffuse light. */
GLfloat light_ambient[] = {0.0, 0.5, 0.0, 1.0};  /* Ambient light. */
GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};  /* Infinite light location. */

bool trackingMouseRotate=false;
bool trackingMouseZoom=false;
bool trackingMouseShift=false;

//Geometry to render
TETRAHEDRON_MESH* tet_mesh=0;
VECTOR<T>* particles;

//Simulation data
BACKWARD_EULER_TIME_STEPPING_3D<T>* backward_euler;
VECTOR<int>* cons_left;
VECTOR<int>* cons_right;
VECTOR<int>* cons;
VECTOR<T>* cons_pos;
bool constrain_both_sides=true;
T boundary_condition_speed=(T)0.1;
bool draw_wire_cube=true;
int wireframe_mode=0;

void processNormalKeys(unsigned char key, int x,int y){
    if(key==27)
        exit(0);
    if(key=='z'||key=='-')
        eye = eye*1.01  ;
    if(key=='Z'||key=='+')
        eye = eye*0.99;
    if(key=='c')
        draw_wire_cube=!draw_wire_cube;
    if(key=='w')
        wireframe_mode=(wireframe_mode+1)%3;
    if(key=='r'){
        tsl=TV();
        eye=def_eye;
        ori=def_ori;
        ACCUMULATED=MATRIX_3X3<T>::Identity();}
}

void processSpecialKeys(int key, int x, int y){
    if(key==GLUT_KEY_UP)   eye = eye*0.99;
    if(key==GLUT_KEY_DOWN) eye = eye*1.01;}

void renderTetMesh(void){
    glClearColor(0.0f,0.0f,0.0f,1.0f); 
    glLoadIdentity();
    gluLookAt(eye(0)+tsl(0),eye(1)+tsl(1),eye(2)+tsl(2),tsl(0),tsl(1),tsl(2),ori(0),ori(1),ori(2));
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 
    glColor3f(0.5,0.5,0.5);
    if(draw_wire_cube){
        glDisable(GL_LIGHT0);
        glDisable(GL_LIGHTING);
        glPushMatrix();
        glTranslatef(tsl(0),0,0);
        glTranslatef(0,tsl(1),0);
        glTranslatef(0,0,tsl(2));
        glutWireCube(1.0);
        glPopMatrix();
        glEnable(GL_LIGHT0);
        glEnable(GL_LIGHTING);}
    glColor3f(1,1,1);
    if(tet_mesh){
        if(tet_mesh->boundary_triangle_mesh){
            for(int tri=0;tri<tet_mesh->boundary_triangle_mesh->Size();tri++){
                int i=(*(tet_mesh->boundary_triangle_mesh))(tri)(0);
                int j=(*(tet_mesh->boundary_triangle_mesh))(tri)(1);
                int k=(*(tet_mesh->boundary_triangle_mesh))(tri)(2);
                TV a((*particles)(3*i),(*particles)(3*i+1),(*particles)(3*i+2));
                TV b((*particles)(3*j),(*particles)(3*j+1),(*particles)(3*j+2));
                TV c((*particles)(3*k),(*particles)(3*k+1),(*particles)(3*k+2));
                TV n=TV::Cross_Product((b-a),(c-a));
                n.Normalize();
                if(wireframe_mode==0 || wireframe_mode==1){
                    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
                    glBegin(GL_TRIANGLES);
                    glNormal3f(n.x(),n.y(),n.z());
                    glVertex3f(a.x(),a.y(),a.z());
                    glVertex3f(b.x(),b.y(),b.z());
                    glVertex3f(c.x(),c.y(),c.z());
                    glEnd();}
                if(wireframe_mode==1 || wireframe_mode==2){
                    glDisable(GL_LIGHT0);
                    glDisable(GL_LIGHTING);
                    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
                    glBegin(GL_TRIANGLES);
                    glVertex3f(a.x(),a.y(),a.z());
                    glVertex3f(b.x(),b.y(),b.z());
                    glVertex3f(c.x(),c.y(),c.z());
                    glEnd();
                    glEnable(GL_LIGHT0);
                    glEnable(GL_LIGHTING);}}}
        glutSwapBuffers();}
}

void track_eye(TV local_u, TV local_v1){
    TV local_v2,local_v3;
    local_v1=ACCUMULATED*local_v1;
    local_v1.Normalize();
    local_u=ACCUMULATED*local_u;
    local_u.Normalize();
    local_v3=TV::Cross_Product(local_v1,local_u);
    local_v3.Normalize();
    local_v2=local_u-local_v1.Dot_Product(local_v1,local_u)*local_v1;
    local_v2.Normalize();
    
    costheta = local_v1.Dot_Product(local_v1, local_u);   //use dot product, local_v1 and local_u have already been normalized
    sintheta = sqrt(1-costheta*costheta);   
    
    MATRIX_3X3<T> Q(local_v1,local_v2,local_v3);
    MATRIX_3X3<T> R(costheta,-sintheta,0,sintheta,costheta,0,0,0,1);
    eye=Q*R*Q.Transposed()*eye;
    ori=Q*R*Q.Transposed()*ori;
    ACCUMULATED=Q*R*Q.Transposed()*ACCUMULATED;}

void startMotionRotate(int x,int y){
    trackingMouseRotate=true;
    v1(0)=rotation_scale*(((T)2*(T)x-(T)winWidth)/(T)winWidth)*one_over_sqrt_3;
    v1(1)=rotation_scale*(((T)winHeight-(T)2*(T)y)/(T)winHeight)*one_over_sqrt_3;}

void stopMotionRotate(int x,int y){
    trackingMouseRotate=false;
    v1(0)=u(0);
    v1(1)=u(1);}

void startMotionZoom(int x,int y){
    trackingMouseZoom=true;
    v1(0)=rotation_scale*(((T)2*(T)x-(T)winWidth)/(T)winWidth)*one_over_sqrt_3;
    v1(1)=rotation_scale*(((T)winHeight-(T)2*(T)y)/(T)winHeight)*one_over_sqrt_3;}

void stopMotionZoom(int x,int y){
    trackingMouseZoom=false;
    v1(0)=u(0);
    v1(1)=u(1);}

void startMotionShift(int x,int y){
    trackingMouseShift=true;
    v1(0)=rotation_scale*(((T)2*(T)x-(T)winWidth)/(T)winWidth)*one_over_sqrt_3;
    v1(1)=rotation_scale*(((T)winHeight-(T)2*(T)y)/(T)winHeight)*one_over_sqrt_3;}

void stopMotionShift(int x,int y){
    trackingMouseShift=false;
    v1(0)=u(0);
    v1(1)=u(1);}

void mouseMotion(int x,int y){
    u(0)=rotation_scale*(((T)2*(T)x-(T)winWidth)/(T)winWidth)*one_over_sqrt_3;
    u(1)=rotation_scale*(((T)winHeight-(T)2*(T)y)/(T)winHeight)*one_over_sqrt_3;

    if(trackingMouseRotate){
        if(fabs(v1(0)-u(0))+fabs(v1(1)-u(1))>0.05) {
            track_eye(u,v1);
            v1(0)=u(0);
            v1(1)=u(1);}} 
    if(trackingMouseZoom){
        if(fabs(v1(1)-u(1))>0.05) {
            eye=eye*exp((u(1)-v1(1))*zoom_scale);
            v1(0)=u(0);
            v1(1)=u(1);}} 
    if(trackingMouseShift){
        if(fabs(v1(0)-u(0))+fabs(v1(1)-u(1))>0.05) {
            tsl=tsl+ACCUMULATED*(v1-u)*eye.Magnitude()*shift_scale;
            v1(0)=u(0);
            v1(1)=u(1);}} 
}

void mouseButton(int button, int state, int x, int y){
    if(button==GLUT_LEFT_BUTTON){
    switch(state){
        case GLUT_DOWN:
            startMotionRotate(x,y);
            break;
        case GLUT_UP:
            stopMotionRotate(x,y);
            break;}}
    if(button==GLUT_RIGHT_BUTTON){
    switch(state){
        case GLUT_DOWN:
            startMotionZoom(x,y);
            break;
        case GLUT_UP:
            stopMotionZoom(x,y);
            break;}}
    if(button==GLUT_MIDDLE_BUTTON){
    switch(state){
        case GLUT_DOWN:
            startMotionShift(x,y);
            break;
        case GLUT_UP:
            stopMotionShift(x,y);
            break;}}
}

void Update_Constrained_Nodes(float time,VECTOR<int>& constrained_left,VECTOR<int>& constrained_right,VECTOR<T>& constrained_locations,const bool left_and_right=true,T speed=(T).25){
    float z_right=1;
    // z_right+=.25*sin(2*PI*time*speed);
    z_right+=.25*(2*PI*time*speed);
    
    float z_left=0;
    // z_left-=.25*sin(2*PI*time*speed);
    z_left-=.25*(2*PI*time*speed);
    
    for(int p=0;p<constrained_right.Size()/3;p++){
        constrained_locations(constrained_left.Size()+3*p+2)=z_right;}
    
    if(left_and_right){
        for(int p=0;p<constrained_left.Size()/3;p++){
            constrained_locations(3*p+2)=z_left;}}}

void Advance_One_Time_Step(){
    Update_Constrained_Nodes(backward_euler->Time(),*cons_left,*cons_right,*cons_pos,constrain_both_sides,boundary_condition_speed);
    backward_euler->Set_Boundary_Conditions(*cons,*cons_pos);
    backward_euler->Advance_One_Time_Step();}

void display(void) {
    Advance_One_Time_Step();
    renderTetMesh();}

void reshape(int width, int height) {
    winWidth=glutGet(GLUT_WINDOW_WIDTH);
    winHeight=glutGet(GLUT_WINDOW_HEIGHT);
    glViewport(0,0,(GLsizei)width,(GLsizei)height);
    glMatrixMode(GL_PROJECTION); 
    glLoadIdentity(); 
    gluPerspective(60,(GLfloat)width/(GLfloat)height,1.0,100.0);
    glMatrixMode(GL_MODELVIEW);}

VECTOR_2D<int> Count_Constrained_Nodes(VECTOR<T>& positions,const bool left_and_right=true){
    float max_z=-FLT_MAX;
    float min_z=FLT_MAX;
    
    for(int p=0;p<positions.Size()/3;p++){
        if(positions(3*p+2)>max_z) max_z=positions(3*p+2);
        if(positions(3*p+2)<min_z) min_z=positions(3*p+2);}
    
    float epsilon=1e-4;
    
    int num_max=0;
    int num_min=0;
        
    for(int p=0;p<positions.Size()/3;p++){
        if(max_z-positions(3*p+2)<epsilon) num_max++;
        if(positions(3*p+2)-min_z<epsilon) num_min++;}
    
    int num_right;
    if(left_and_right) num_right=3*num_max;
    else num_right=0;
    
    return VECTOR_2D<int>(3*num_min,num_right);}

void Initialize_Constrained_Nodes(VECTOR<int>& constrained_left,VECTOR<int>& constrained_right,VECTOR<T>& positions,const bool left_and_right=true){
    float max_z=-FLT_MAX;
    float min_z=FLT_MAX;
    
    for(int p=0;p<positions.Size()/3;p++){
        if(positions(3*p+2)>max_z) max_z=positions(3*p+2);
        if(positions(3*p+2)<min_z) min_z=positions(3*p+2);}
    
    float epsilon=1e-4;
    
    int num_max=0;
    int num_min=0;
    
    for(int p=0;p<positions.Size()/3;p++){
        if(max_z-positions(3*p+2)<epsilon && left_and_right){constrained_right(num_max++)=3*p;constrained_right(num_max++)=3*p+1;constrained_right(num_max++)=3*p+2;}
        if(positions(3*p+2)-min_z<epsilon){constrained_left(num_min++)=3*p;constrained_left(num_min++)=3*p+1;constrained_left(num_min++)=3*p+2;}}
    
    assert(num_max==constrained_right.Size());
    assert(num_min==constrained_left.Size());}

int main(int argc, char **argv)
{
    // SET UP GEOMETRY
    const int n=7; //the number of particles per side, n>=2 
    
    // NUMBER THE PARTICLES
    T d_x=(T)1.0/(T)(n-1);   
    GRID_3D<T> grid(n,d_x,(T)0,(T)0,(T)0);
    DEFORMABLE_OBJECT_3D<T> deformable_object(grid);
    tet_mesh=&(deformable_object.Tetrahedron_Mesh());
    tet_mesh->Initialize_Oriented_Boundary_Triangles();
    particles=&(deformable_object.Positions());
    for(int k=0;k<n;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<n;i++){
                (*particles)(3*(i+n*j+n*n*k))=i*d_x;
                (*particles)(3*(i+n*j+n*n*k)+1)=j*d_x;
                (*particles)(3*(i+n*j+n*n*k)+2)=k*d_x;}}}
    
    // ELASTIC CONSTITUTIVE MODEL
    // LINEAR_ELASTICITY_3D<T> le((float)10000,(float).3);
    FIXED_COROTATED_3D<T> le((float)10000,(float).3);

    FEM_HYPERELASTICITY_3D<T> fem(deformable_object.Tetrahedron_Mesh(),deformable_object.Positions());
    fem.Set_Constitutive_Model(le);
    fem.Initialize_Undeformed_Configuration();
    
    // CHOOSE TIME STEPPING: BACKWARD EULER
    T dt=(T)1/(T)300;
    int number_steps=30;
    float start_time=(float)0;
    float end_time=(float)number_steps*dt;
    
    BACKWARD_EULER_TIME_STEPPING_3D<T> be(dt,end_time,start_time,deformable_object);
    backward_euler=&be;
    be.Set_Elastic_Forces(fem);
    VECTOR<T> nodal_volumes(deformable_object.Positions().Size());
    fem.Nodal_Volume_Fractions(nodal_volumes);
    
    be.Initialize_BE_Matrix(nodal_volumes);
    be.Initialize_CG();
    
    // BOUNDARY CONDITIONS
    VECTOR_2D<int> numbers_constrained=Count_Constrained_Nodes(*particles,constrain_both_sides);
    VECTOR<int> constrained_left(numbers_constrained.x());
    cons_left=&constrained_left;
    VECTOR<int> constrained_right(numbers_constrained.y());
    cons_right=&constrained_right;
    Initialize_Constrained_Nodes(constrained_left,constrained_right,*particles,constrain_both_sides);       
    
    int number_constrained=constrained_left.Size()+constrained_right.Size();
    VECTOR<int> constrained(number_constrained);
    cons=&constrained;
    for(int i=0;i<constrained_left.Size();i++) constrained(i)=constrained_left(i);
    for(int i=0;i<constrained_right.Size();i++) constrained(i+constrained_left.Size())=constrained_right(i);
    VECTOR<T> constrained_locations(number_constrained);
    for(int i=0;i<number_constrained;i++) constrained_locations(i)=(*particles)(constrained(i));
    cons_pos=&constrained_locations;
    
    glutInit(&argc,argv); 
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_DEPTH); 
    glutInitWindowSize(1024,768); 
    glutInitWindowPosition(100,100);
    glutCreateWindow("Simulation Visualization"); 
    
    winWidth=glutGet(GLUT_WINDOW_WIDTH);
    winHeight=glutGet(GLUT_WINDOW_HEIGHT);
    
    glutDisplayFunc(display); 
    glutIdleFunc(display); 
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(processNormalKeys);
    glutSpecialFunc(processSpecialKeys);
    glutReshapeFunc(reshape);
    
    glEnable(GL_DEPTH_TEST); 

    glLightfv(GL_LIGHT0,GL_DIFFUSE,light_diffuse);
    glLightfv(GL_LIGHT0,GL_AMBIENT,light_ambient);
    glLightfv(GL_LIGHT0,GL_POSITION,light_position);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glPolygonOffset(1,1);
    glEnable(GL_POLYGON_OFFSET_FILL);
   
    glutMainLoop(); 
    
    if(tet_mesh) delete tet_mesh;
    if(particles) delete particles;
    return 0;
}


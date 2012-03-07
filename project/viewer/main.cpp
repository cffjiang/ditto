#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>
#include <GL/glut.h>
#include <ditto/public_library/visualization/arcball.h>
#include <ditto/public_library/geometry/simplex.h>
#include <ditto/public_library/algebra/linear_algebra.h>

using namespace ditto::visualization;

typedef double T;
std::string prefix;
int frame = 1;
bool auto_update = false;

typedef ditto::algebra::VECTOR_3D<T> node_type;
typedef std::vector<node_type> node_list_type;
typedef ditto::algebra::VECTOR_3D<int> element_type;
typedef std::vector<element_type> element_list_type;

struct Mesh {
    node_list_type nodes;
    element_list_type elements;
    node_list_type vertex_normals;
} mesh;

struct Ball {
    T x;
    T y;
    T z;
    T r;
};

std::vector<Ball> ball_list;

GLint uAmbient, uDiffuse, uSpecular, uLightPos, uShininess;
static float aspect_ratio = 1.0f;
static int width, height;
vec eye( 0.0f, 0.0f, -20.0f );
vec centre( 0.0f, 0.0f, 0.0f );
vec object_translate( 0.0f, 0.0f, 0.0f );
const vec up( 0.0f, 1.0f, 0.0f );
const float SPHERE_RADIUS = 5.0f;
const int SPHERE_LAT_SLICES = 12;
const int SPHERE_LONG_SLICES = 24;
const float PI = 3.141592654f;

void read_new_frame(); 
static void reset_view(int w, int h);
inline float randf();
static void shutdown_scene();
inline vec rotate_x( vec v, float sin_ang, float cos_ang );
inline vec rotate_y( vec v, float sin_ang, float cos_ang );
static void draw_sphere();
static void click_scene(int x, int y);
static void drag_scene(int x, int y);
static void draw_scene();
static void resize(int w, int h);
static void display();
static void key(unsigned char key, int x, int y);
static void idle();
static void mouse_button(int button, int state, int x, int y);
static void mouse_motion(int x, int y);

void read_new_frame() {
    std::stringstream ssframe;
    ssframe << frame;
    std::string filename = prefix + ssframe.str();

    std::ifstream in(filename.c_str());

    if (in == NULL) return;

    std::cout << "Reading Frame " << frame << std::endl;

    std::istringstream ss;
    std::string s;

    int node_count = 0;
    int face_count = 0;
    int num_nodes = -1;
    int num_faces = -1;
    bool reading_nodes = false;
    bool reading_faces = false;  

    while (getline(in, s)) {
        ss.clear();
        ss.str(s);
        if (reading_nodes == true) {
            ss >> mesh.nodes[node_count](0) >> mesh.nodes[node_count](1) >> mesh.nodes[node_count](2);
            node_count++;
            if (node_count == num_nodes) reading_nodes = false;
        }
        
        else if (reading_faces == true) {
            ss >> mesh.elements[face_count](0) >> mesh.elements[face_count](1) >> mesh.elements[face_count](2);
            face_count++;
            if (face_count == num_faces) reading_faces = false;       
        }
        
        else if (s.substr(0,5) == "NODES") {
            std::string tmp1;
            ss >> tmp1 >> num_nodes;
            if (mesh.nodes.size() != num_nodes) {
                mesh.nodes.resize(num_nodes); }
            reading_nodes = true;
            reading_faces = false;
        }
        
        else if (s.substr(0,4) == "TRIS") {
            std::string tmp1, tmp2;
            ss >> tmp1 >> num_faces;
            if (mesh.elements.size() != num_faces) {
                mesh.elements.resize(num_faces); }
            reading_nodes = false;            
            reading_faces = true;            
        }    

        else if (s.substr(0,4) == "BALL") {
            ball_list.clear();
            std::string tmp1;
            Ball ball;
            ss >> tmp1 >> ball.x >> ball.y >> ball.z >> ball.r;
            ball_list.push_back(ball);
         }   
    }

    // compute normals
    mesh.vertex_normals.clear();
    mesh.vertex_normals.resize(mesh.nodes.size());
    for (unsigned int i=0; i<mesh.elements.size(); i++) {
        int node1 = mesh.elements[i](0);
        int node2 = mesh.elements[i](1);
        int node3 = mesh.elements[i](2);

        node_type edge1 = mesh.nodes[node2] - mesh.nodes[node1];
        node_type edge2 = mesh.nodes[node3] - mesh.nodes[node1];
        node_type pre_normal = edge1.Cross_Product(edge2);
        node_type normal = pre_normal  * (1.0 / pre_normal.Magnitude());

        mesh.vertex_normals[node1] =  mesh.vertex_normals[node1] + normal;
        mesh.vertex_normals[node2] =  mesh.vertex_normals[node2] + normal;
        mesh.vertex_normals[node3] =  mesh.vertex_normals[node3] + normal;
    }

     for (unsigned int i=0; i<mesh.vertex_normals.size(); i++) {
         T mag = mesh.vertex_normals[i].Magnitude();
         mesh.vertex_normals[i] =  mesh.vertex_normals[i] * (1.0/mag);
     }


}

static void reset_view(int w, int h){
    width = w;
    height = h;
    aspect_ratio = (float) width / (float) height;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective( 50.0f, aspect_ratio, 1.0f, 50.0f );
    gluLookAt(
        eye.x, eye.y, eye.z,
        centre.x, centre.y, centre.z,
        up.x, up.y, up.z );
    arcball_setzoom( SPHERE_RADIUS, eye, up );
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ; }

void set_colour(float r, float g, float b){
    float ambient = 0.2f;
    float diffuse = 0.7f;
    float specular = 0.7f;
    GLfloat mat[4];
    mat[0] = ambient*r;
    mat[1] = ambient*g;
    mat[2] = ambient*b;
    mat[3] = 1.0;
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, mat);
    mat[0] = diffuse*r;
    mat[1] = diffuse*g;
    mat[2] = diffuse*b;
    mat[3] = 1.0;
    glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, mat);
    mat[0] = specular*r;
    mat[1] = specular*g;
    mat[2] = specular*b;
    mat[3] = 1.0;
    glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat);
    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 1.0);}

void my_init(){
    GLfloat ambient[] = { 0.1, 0.1, 0.1, 1.0 };

    GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat position[] = { 10.0, 0.0, -30.0, 1.0 };

    GLfloat diffuse2[] = { 0.3, 0.3, 0.3, 1.0 };
    GLfloat specular2[] = { 0.3, 0.3, 0.3, 1.0 };
    GLfloat position2[] = { 0.0, 100.0, 00.0, 1.0 };

    GLfloat diffuse3[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat specular3[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat position3[] = { 10.0, 0.0, 30.0, 1.0 };

    GLfloat lmodel_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat local_view[] = { 0.0 };
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    glLightf(GL_LIGHT0, GL_SHININESS, 100) ;
    glLightfv(GL_LIGHT0, GL_POSITION, position);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE) ;
    glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse2);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular2);
    glLightfv(GL_LIGHT1, GL_POSITION, position2);
    glLightf(GL_LIGHT1, GL_SHININESS, 500) ;

    glLightfv(GL_LIGHT2, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuse3);
    glLightfv(GL_LIGHT2, GL_SPECULAR, specular3);
    glLightf(GL_LIGHT2, GL_SHININESS, 100) ;
    glLightfv(GL_LIGHT2, GL_POSITION, position3);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glPixelStorei(GL_PACK_ALIGNMENT,1) ;
    glPixelStorei(GL_UNPACK_ALIGNMENT,1) ;
    glShadeModel(GL_SMOOTH) ;

    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);  
    
    read_new_frame();}

inline float randf(){
    return ((1.0f / 127.f) * (((float)(rand() & 255)) - 128.0f)) ;}

static void shutdown_scene(){}

inline vec rotate_x( vec v, float sin_ang, float cos_ang ){
    return vec(
        v.x,
        (v.y * cos_ang) + (v.z * sin_ang),
        (v.z * cos_ang) - (v.y * sin_ang)
    );}

inline vec rotate_y( vec v, float sin_ang, float cos_ang ){
    return vec(
        (v.x * cos_ang) + (v.z * sin_ang),
        v.y,
	    (v.z * cos_ang) - (v.x * sin_ang)
    );}

static void draw_sphere() {
    glutSolidSphere(1,30,20);}

static void click_scene(int x, int y){
    int invert_y = (height - y) - 1;
    arcball_start(x,invert_y);}

static void drag_scene(int x, int y){
    int invert_y = (height - y) - 1;
    arcball_move(x,invert_y);}

static void draw_scene(){
    glTranslatef( object_translate.x, object_translate.y, object_translate.z );
    arcball_rotate();

    // draw things below
    set_colour( 0.0f, 1.0f, 0.0f );

    // draw mesh

    for (int i=0; i<mesh.elements.size(); i++) {
        glBegin(GL_TRIANGLES);

        glNormal3f( mesh.vertex_normals[mesh.elements[i](0)](0),  mesh.vertex_normals[mesh.elements[i](0)](1),  mesh.vertex_normals[mesh.elements[i](0)](2));
        glVertex3f( mesh.nodes[mesh.elements[i](0)](0),  mesh.nodes[mesh.elements[i](0)](1),  mesh.nodes[mesh.elements[i](0)](2));

        glNormal3f( mesh.vertex_normals[mesh.elements[i](1)](0),  mesh.vertex_normals[mesh.elements[i](1)](1),  mesh.vertex_normals[mesh.elements[i](1)](2));
        glVertex3f( mesh.nodes[mesh.elements[i](1)](0),  mesh.nodes[mesh.elements[i](1)](1),  mesh.nodes[mesh.elements[i](1)](2));

        glNormal3f( mesh.vertex_normals[mesh.elements[i](2)](0),  mesh.vertex_normals[mesh.elements[i](2)](1),  mesh.vertex_normals[mesh.elements[i](2)](2));
        glVertex3f( mesh.nodes[mesh.elements[i](2)](0),  mesh.nodes[mesh.elements[i](2)](1),  mesh.nodes[mesh.elements[i](2)](2));

        glEnd();
    }


    // draw ball
    if (ball_list.size()!=0) {
        set_colour(1.0, 0.0, 0.0);
        for (unsigned int i=0; i<ball_list.size(); i++) {
            glPushMatrix();
            glTranslatef(ball_list[i].x, ball_list[i].y, ball_list[i].z);
            glutSolidSphere(ball_list[i].r,30,20);
            glPopMatrix();} }
}

static void resize(int w, int h){
    reset_view(w,h);}

static void display(){
    if (auto_update == true) {
        frame++;
        read_new_frame();}

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    draw_scene();
    glutSwapBuffers();}

static void key(unsigned char key, int x, int y){
    switch (key){
        case 27 : 
        case 'q':
            shutdown_scene();
            exit(0);
            break;
        case 'w':
            object_translate.y += .2;
            break;
        case 's':
            object_translate.y -= .2;
            break;
        case 'a':
            object_translate.x += .2;
            break;
        case 'd':
            object_translate.x -= .2;
            break;
        case '-':
            object_translate.z += 1;
            break;
        case '+':
            object_translate.z -= 1;
            break;
        case 'n':
            frame++;
            read_new_frame();
            break;
        case 'p':
            auto_update = !auto_update;
            break;
        case 'r':
            frame = 1;
            read_new_frame();
            break;
        case '1':
            glEnable(GL_LIGHTING);
            glPolygonMode(GL_FRONT, GL_FILL);
            glPolygonMode(GL_BACK, GL_FILL);  
            break;
        case '2':
            glDisable(GL_LIGHTING);
            glColor3f(0.0, 1.0, 0.0);
            glPolygonMode(GL_FRONT, GL_LINE);
            glPolygonMode(GL_BACK, GL_LINE);  
            break;
        default :
            break;}
    glutPostRedisplay();}

static void idle(){
    glutPostRedisplay();}

static void mouse_button(int button, int state, int x, int y){
    if ( state == GLUT_DOWN ) click_scene(x,y);}

static void mouse_motion(int x, int y){
    drag_scene(x,y);}

int main(int argc, char ** argv){
    prefix = argv[1];
    arcball_reset();
    glutInit(&argc, argv);
    glutInitWindowSize(800,600);
    glutInitWindowPosition(10,10);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Ditto OpenGL Viewer");
    my_init();
    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutIdleFunc(idle);
    glutMouseFunc(mouse_button);
    glutMotionFunc(mouse_motion);
    glClearColor(0.0f,0.0f,0.0f,1.0f);
    //glEnable(GL_CULL_FACE);
    //glCullFace(GL_BACK);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glutMainLoop();
    return 0;}

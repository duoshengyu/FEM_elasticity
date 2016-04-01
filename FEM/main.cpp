#include "App.h"
#include "camera.h"
#include <vector>
#include <map>
#include "GameTimer.h"
#include <GL/wglew.h>
#include <GL/freeglut.h>
#include "deformableobject.h"
#include "hash.h"

#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "OpenMeshCore.lib")
#pragma comment(lib, "OpenMeshTools.lib")
//#pragma pack(push,1)

//-----------------------------------------------------------------------
//two debug file 

extern ofstream debug;
extern ofstream debug2;
float scale = 3.0f;

const int width = 1024, height = 768;

float timeStep = 1 / 60.0f;
float currentTime = 0;
float accumulator = timeStep;

int isMouseButtonDown = 0;


int oldX = 0, oldY = 0;
float rX = 15, rY = 0;
int state = 1;
float dist = -10.0f;
const int GRID_SIZE = 5;


GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

glm::vec3 Up = glm::vec3(0, 1, 0), Right, viewDir;

LARGE_INTEGER frequency;        // ticks per second
LARGE_INTEGER t1, t2;           // ticks
float frameTimeQP = 0;
float frameTime = 0;

float startTime = 0, fps = 0;
int totalFrames = 0;

char info[MAX_PATH] = { 0 };


void DrawGrid()
{
	glBegin(GL_QUADS);
	glColor3f(0.5f, 0.5f, 0.5f);
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
	{
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(-5.0f + j, 0, -5.0f + i);
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(-5.0f + j + 1, 0, -5.0f + i);
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(-5.0f + j + 1, 0, -5.0f + i + 1);
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(-5.0f + j, 0, -5.0f + i + 1);
	}
	glEnd();
}

//large prime number use for hashing 
__int64 p1 = 73856093;
__int64 p2 = 19349663;
__int64 p3 = 83492791;
float g_l = 0;//global grid size
int n = 199;//hash table size
float a = 9.81;//penalty coefficient
HashMap g_HashTable(n);
vector<DeformableObject> g_models;

vector<vector<F_add>> g_F;



namespace shader
{
	GLuint loadShader(const char * filename, GLenum shader_type, bool check_errors)
	{
		GLuint result = 0;
		string shaderCode = "";
		ifstream codeStream(filename, ios::in);

		if (codeStream.is_open())
		{
			string Line = "";
			while (getline(codeStream, Line))
			{
				shaderCode += Line;
				shaderCode += "\n";
			}

			codeStream.close();
		}
		else
		{
			cout << filename << " can not open! " << endl;
		}
		result = glCreateShader(shader_type);

		// Compile Vertex Shader
		//cout << "Compiling shader : " << filename << endl;
		char const * codePointer = shaderCode.c_str();

		glShaderSource(result, 1, &codePointer, NULL);
		glCompileShader(result);

		if (check_errors)
		{
			GLint status = 0;
			glGetShaderiv(result, GL_COMPILE_STATUS, &status);

			if (!status)
			{
				char buffer[4096];
				glGetShaderInfoLog(result, 4096, NULL, buffer);

				cout << filename << buffer << endl;
			}
		}

		return result;
	}

}


class FEMTest : public App
{
public:
	FEMTest();
	~FEMTest();

	bool                    Init();
	void                    UpdateScene();
	void                    Rendering();
	void                    onResize(GLFWwindow* window, int w, int h);

	void                    onMouseWheel(GLFWwindow* window, double x, double y);
	void                    onMouseMove(GLFWwindow* window, double xd, double yd);
	void                    onMouseButton(GLFWwindow* window, int button, int action, int mods);
	void                    onKey(GLFWwindow* window, int key, int scancode, int action, int mods);

private:
	void                    buildGeometryBuffers();
	void                    buildShader();
	void                    setGridCellSize();

private:

	GLuint                  vao;
	GLuint                  program;
	GLuint                  mvp_matrix;
	GLuint                  m_matrix;
	GLuint                  v_matrix;
	GLuint                  l_position;
	TrackballCamera         camera;

	GameTimer                timer;
	DeformableObject         bunny;
	DeformableObject         bunny2;
	DeformableObject         bunny3;
	//DeformableObject         block;
	HashMap                  HashTable;
};

int main(void)
{
	FEMTest *theApp = new FEMTest;

	if (!theApp->Init())
		return 0;
	theApp->Run();
	delete theApp;
	return 0;
}

FEMTest::FEMTest() : App(), timer(), HashTable(n), camera(width, height),
bunny(0.1f, 0.30f, 0.0f, false), bunny2(0.0f, 0.0f, 0.0f, true), bunny3(-0.03f, 0.4f, 0.0f, false)//, block(10, 3, 3, 0.1f, 0.1f, 0.1f),
{
	//g_models.push_back(block);
	g_models.push_back(bunny);
	g_models.push_back(bunny2);
	g_models.push_back(bunny3);
	g_F.resize(g_models.size());
	mWidth = width;
	mHeight = height;
}

FEMTest::~FEMTest()
{
	// Cleanup VBO and shader
}

bool FEMTest::Init()
{
	if (!App::Init())
		return false;

	//startTime = (float)glutGet(GLUT_ELAPSED_TIME);
	//currentTime = startTime;

	// get ticks per second
	//QueryPerformanceFrequency(&frequency);

	// start timer
	//QueryPerformanceCounter(&t1);

	timer.Reset();
	startTime = timer.getCurrenTime();
	currentTime = startTime;

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	wglSwapIntervalEXT(0);

	onResize(window, width, height);

	setGridCellSize();

	buildShader();

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	return true;
}

void FEMTest::onResize(GLFWwindow* window, int nw, int nh)
{
	App::onResize(window, nw, nh);

	glViewport(0, 0, nw, nh);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//gluPerspective(60, (GLfloat)nw / (GLfloat)nh, 0.1f, 100.0f);
	//
	//glGetIntegerv(GL_VIEWPORT, viewport);
	//glGetDoublev(GL_PROJECTION_MATRIX, P);
	//
	//glMatrixMode(GL_MODELVIEW);
	camera.SetProj(60, (GLfloat)nw / (GLfloat)nh, 0.1f, 100.0f);
	camera.SetView(glm::vec3(0.0f, 1.0f, 2.0f), glm::vec3(0.0f, 0.0f, 0.0f));
}
//-----------------------------------------------------------------------
//Do SP collision search
//-----------------------------------------------------------------------
void FEMTest::UpdateScene()
{
	if (accumulator >= timeStep)
	{
		HashTable.T = timer.TotalTime();

		for (int i = 0; i < g_models.size(); ++i)
		{
			g_models[i].firstPass(&HashTable, i);
		}
		for (int i = 0; i < g_F.size(); ++i)
		{
			g_F[i].clear();
		}
		for (int i = 0; i < g_models.size(); ++i)
		{
			for (int j = 0; j < g_models[i].total_points; ++j)
			{
				g_models[i].F[j].x = 0;
				g_models[i].F[j].y = 0;
				g_models[i].F[j].z = 0;
			}
		}
		for (int i = 0; i < g_models.size(); ++i)
		{
			g_models[i].secondPass(&HashTable, i);
		}

		for (int i = 0; i < g_models.size(); ++i)
		{
			g_models[i].F_ADD(i);
		}

		for (int i = 0; i < g_models.size(); ++i)
		{
            g_models[i].StepPhysics(1/500.0f);
		}
		for (int i = 0; i < g_models.size(); ++i)
		{
			g_models[i].CalculateVN();
		}
			
		accumulator -= timeStep;
	}
}
//-----------------------------------------------------------------------
//draw
//-----------------------------------------------------------------------
void FEMTest::Rendering()
{
	timer.Tick();
	float newTime = timer.getCurrenTime();
	frameTime = newTime - currentTime;
	currentTime = newTime;
	//accumulator += frameTime;

	//Using high res. counter
	//QueryPerformanceCounter(&t2);
	// compute and print the elapsed time in millisec
	//frameTimeQP = (t2.QuadPart - t1.QuadPart) * 1000.0f / frequency.QuadPart;
	frameTimeQP = timer.DeltaTime();
	//t1 = t2;
	accumulator += frameTimeQP;

	++totalFrames;
	if ((newTime - startTime) > 1000)
	{
		float elapsedTime = (newTime - startTime);
		fps = (totalFrames / elapsedTime) * 1000;
		startTime = newTime;
		totalFrames = 0;
	}

	sprintf_s(info, "FPS: %3.2f, Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f, Stiffness Warp: %s", fps, frameTime, frameTimeQP, bunny.bUseStiffnessWarping ? "On" : "Off");
	glfwSetWindowTitle(window, info);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);    
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//draw grid
	DrawGrid();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glm::mat4 MVP = camera.getMVP();
    glm::mat4 M = camera.getM();
    glm::mat4 V = camera.getV();
    
    glUseProgram(program);
    
    glUniformMatrix4fv(mvp_matrix, 1, GL_FALSE, glm::value_ptr(MVP));
    glUniformMatrix4fv(v_matrix, 1, GL_FALSE, glm::value_ptr(V));
    glUniformMatrix4fv(m_matrix, 1, GL_FALSE, glm::value_ptr(M));
    glm::vec3 lightPos = glm::vec3(0, 0.25, 0.5);
    glUniform3f(l_position, lightPos.x, lightPos.y, lightPos.z);

	//glMatrixMode(GL_MODELVIEW);
	//glScaled(2.0f, 2.0f, 2.0f);

		for (int i = 0; i < g_models.size(); ++i)
			g_models[i].RenderModel();
	/*
	//draw grid
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (int i = 0; i < g_models[0].gi - 1; ++i)
	//for (int i = 2; i < 3; ++i)
		//for (int j = 7; j < 8; ++j)
		for (int j = 0; j < g_models[0].gj - 1; ++j)
		//for (int k = 1; k < 2; ++k)
			for (int k = 0; k < g_models[0].gk - 1; ++k)
			{
		int n1 = i * g_models[0].gj * g_models[0].gk + j * g_models[0].gk + k;
		int n2 = n1 + 1;
		int n3 = i * g_models[0].gj * g_models[0].gk + (j + 1) * g_models[0].gk + k;
		int n4 = n3 + 1;
		int n5 = (i + 1) * g_models[0].gj * g_models[0].gk + j * g_models[0].gk + k;
		int n6 = n5 + 1;
		int n7 = (i + 1) * g_models[0].gj * g_models[0].gk + (j + 1) * g_models[0].gk + k;
		int n8 = n7 + 1;



		glm::vec3 p1 = g_models[0].grid[n1];
		glm::vec3 p2 = g_models[0].grid[n2];
		glm::vec3 p3 = g_models[0].grid[n3];
		glm::vec3 p4 = g_models[0].grid[n4];
		glm::vec3 p5 = g_models[0].grid[n5];
		glm::vec3 p6 = g_models[0].grid[n6];
		glm::vec3 p7 = g_models[0].grid[n7];
		glm::vec3 p8 = g_models[0].grid[n8];

		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p3.x, p3.y, p3.z);
		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p3.x, p3.y, p3.z);

		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p5.x, p5.y, p5.z);
		glVertex3f(p2.x, p2.y, p2.z);		glVertex3f(p6.x, p6.y, p6.z);
		glVertex3f(p7.x, p7.y, p7.z);		glVertex3f(p3.x, p3.y, p3.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p8.x, p8.y, p8.z);

		glVertex3f(p5.x, p5.y, p5.z);		glVertex3f(p6.x, p6.y, p6.z);
		glVertex3f(p5.x, p5.y, p5.z);		glVertex3f(p7.x, p7.y, p7.z);
		glVertex3f(p6.x, p6.y, p6.z);		glVertex3f(p8.x, p8.y, p8.z);
		glVertex3f(p7.x, p7.y, p7.z);		glVertex3f(p8.x, p8.y, p8.z);
			}
	glEnd();
	
	
	{
	//draw gradient of grid
	glColor3f(0.0, 0.5, 1.0);
	glBegin(GL_LINES);
	for (int i = 1; i < g_models[0].gi - 1; ++i)
	//for (int i = 2; i < 4; ++i)
		//for (int j = 7; j < 9; ++j)
	for (int j = 1; j < g_models[0].gj - 1; ++j)
	//for (int k = 1; k < 3; ++k)
		for (int k = 1; k < g_models[0].gk - 1; ++k)
		{
		int index = i * g_models[0].gj * g_models[0].gk + j * g_models[0].gk + k;
		glm::vec3 p1, p2;
		p1 = g_models[0].grid[index];
		p2 = g_models[0].gradientFild[index];
		//p2 = glm::normalize(p2);
		p2 = p1 + p2;
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
		}
	glEnd();
}
	*/
	/*
	//draw gradient of vertex
	glColor3f(1.0, 1.0, 0.0);
	glBegin(GL_LINES);
	//for (int i = 30; i < 31; ++i)
	for (int i = 0; i < g_models[0].total_points; ++i)
	{
		glm::vec3 p1, p2;
		p1 = g_models[0].X[i];
		p2 = g_models[0].gradientV[i];
		//p2 = glm::normalize(p2);
		p2 = p1 + p2;	
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
	}
	glEnd();*/
	/*
	for (int i = 30; i < 31; ++i)
	{
		for (int numT = 0; numT < g_models[0].tetrahedraOfVertices[i].size(); ++numT)
		{
			int tetrahedraIndex = g_models[0].tetrahedraOfVertices[i][numT];
			glColor3f(0.75, 0.75, 0.75);
			{
				glBegin(GL_LINES);
				int n0 = g_models[0].tetrahedra[tetrahedraIndex].indices[0];
				int n1 = g_models[0].tetrahedra[tetrahedraIndex].indices[1];
				int n2 = g_models[0].tetrahedra[tetrahedraIndex].indices[2];
				int n3 = g_models[0].tetrahedra[tetrahedraIndex].indices[3];
				glm::vec3 p1 = glm::vec3(g_models[0].X[n0].x, g_models[0].X[n0].y, g_models[0].X[n0].z);
				glm::vec3 p2 = glm::vec3(g_models[0].X[n1].x, g_models[0].X[n1].y, g_models[0].X[n1].z);
				glm::vec3 p3 = glm::vec3(g_models[0].X[n2].x, g_models[0].X[n2].y, g_models[0].X[n2].z);
				glm::vec3 p4 = glm::vec3(g_models[0].X[n3].x, g_models[0].X[n3].y, g_models[0].X[n3].z);

				glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p1.x, p1.y, p1.z);
				glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p2.x, p2.y, p2.z);
				glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p3.x, p3.y, p3.z);

				glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p2.x, p2.y, p2.z);
				glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p3.x, p3.y, p3.z);

				glVertex3f(p2.x, p2.y, p2.z);		glVertex3f(p3.x, p3.y, p3.z);
				glEnd();
			}
		}
	}
	*/
	/*
	//is point in tetrahedron

		glColor3f(0.75, 0.75, 0.75);
		{
		glBegin(GL_LINES);

		glm::vec3 p1 = glm::vec3(-0.0158555,   0.149626, - 0.0252782);
		glm::vec3 p2 = glm::vec3(-0.0208207,   0.141843, - 0.0189655);
		glm::vec3 p3 = glm::vec3(-0.0289981,   0.135866, - 0.0183379);
		glm::vec3 p4 = glm::vec3(-0.0131466,   0.140465, - 0.0290179);

			glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p1.x, p1.y, p1.z);
			glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p2.x, p2.y, p2.z);
			glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p3.x, p3.y, p3.z);

			glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p2.x, p2.y, p2.z);
			glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p3.x, p3.y, p3.z);

			glVertex3f(p2.x, p2.y, p2.z);		glVertex3f(p3.x, p3.y, p3.z);
		glEnd();
}
*/
}
void FEMTest::onMouseWheel(GLFWwindow* window, double x, double y)
{

	float scale = 1.0f;
	float mouseWheelScale = 0.1f;
	scale += mouseWheelScale  * (float)y;
	//camera.setScaleFactor(scale);
	//camera.setMmworldScle();
}
void FEMTest::onMouseMove(GLFWwindow* window, double xd, double yd)
{
	double x = xd;
	double y = yd;
	if (camera.IsMouseLButtonDown())
	{
		glfwGetCursorPos(window, &xd, &yd);
		camera.SetCurMousePosition(xd, yd);
		camera.computeQuat();
		camera.setMmworldQuat();
	}
	if (isMouseButtonDown)
	{
		if (bunny.GetSelectIndex() == -1) {
			if (state == 0)
				dist *= (1 + (y - oldY) / 60.0f);
			else
			{
				rY += (x - oldX) / 5.0f;
				rX += (y - oldY) / 5.0f;
			}
		}
		else {

			float delta = 1000 / abs(dist);
			float valX = (x - oldX) / delta;
			float valY = (oldY - y) / delta;

			bunny.V[bunny.GetSelectIndex()] = glm::vec3(0);
			bunny.X[bunny.GetSelectIndex()].x += Right[0] * valX;
			float newValue = bunny.X[bunny.GetSelectIndex()].y + Up[1] * valY;
			if (newValue>0)
				bunny.X[bunny.GetSelectIndex()].y = newValue;
			bunny.X[bunny.GetSelectIndex()].z += Right[2] * valX + Up[2] * valY;
		}
		oldX = x;
		oldY = y;
	}

}
void FEMTest::onMouseButton(GLFWwindow* window, int button, int action, int mods)
{
	double xd, yd;

	if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_PRESS))
	{
		camera.SetMouseLButtonStat(true);
		glfwGetCursorPos(window, &xd, &yd);
		camera.initMousePosition(xd, yd);


		glfwGetCursorPos(window, &xd, &yd);
		oldX = xd;
		oldY = yd;
		int window_y = (height - yd);
		float norm_y = float(window_y) / float(height / 2.0);
		int window_x = xd;
		float norm_x = float(window_x) / float(width / 2.0);

		float winZ = 0;
		glReadPixels(xd, height - yd, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
		double objX = 0, objY = 0, objZ = 0;

		gluUnProject(window_x, window_y, winZ, MV, P, viewport, &objX, &objY, &objZ);

		glm::vec3 pt(objX, objY, objZ);
		printf("\nObj [ %3.3f,%3.3f,%3.3f ]", objX, objY, objZ);
		size_t i = 0;
		for (i = 0; i < bunny.total_points; i++) {
			if (glm::distance(bunny.X[i], pt)<0.01) {
				bunny.SetSelectIndex(i);

				printf("Intersected at %d\n", i);
				printf("Pt [ %3.3f,%3.3f,%3.3f ]\n", bunny.X[i].x, bunny.X[i].y, bunny.X[i].z);
				break;
			}
		}
		isMouseButtonDown = 1;
	}
	else if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_RELEASE))
	{
		camera.SetMouseLButtonStat(false);

		bunny.SetSelectIndex(-1);
		bunny.UpdateOrientation();
		isMouseButtonDown = 0;
	}

	if ((button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_PRESS))
	{
		state = 0;
	}
	else if ((button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_RELEASE))
	{
		state = 1;
	}

}
void FEMTest::onKey(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	App::onKey(window, key, scancode, action, mods);
	if ((key == GLFW_KEY_SPACE) && (action == GLFW_PRESS))
	{
		//bunny.bUseStiffnessWarping = !bunny.bUseStiffnessWarping;
		//printf("Stiffness Warping %s\n", bunny.bUseStiffnessWarping ? "On" : "Off");
		for (int i = 0; i < g_models.size(); ++i)
			g_models[i].Reset();
	}

}
void FEMTest::buildGeometryBuffers()
{

}
//-----------------------------------------------------------------------
//read shader
//-----------------------------------------------------------------------
void FEMTest::buildShader()
{
	GLuint vs;
	GLuint fs;

	vs = shader::loadShader("lighting.vs.glsl", GL_VERTEX_SHADER, true);
	fs = shader::loadShader("lighting.fs.glsl", GL_FRAGMENT_SHADER, true);

	//vs = shader::loadShader("ColorVertexShader.glsl", GL_VERTEX_SHADER, true);
	//fs = shader::loadShader("ColorFragmentShader.glsl", GL_FRAGMENT_SHADER, true);
	if (program)
		glDeleteProgram(program);

	program = glCreateProgram();
	glAttachShader(program, vs);
	glAttachShader(program, fs);

	glLinkProgram(program);

	glDeleteShader(vs);
	glDeleteShader(fs);

	mvp_matrix = glGetUniformLocation(program, "MVP");
	v_matrix = glGetUniformLocation(program, "V");
	m_matrix = glGetUniformLocation(program, "M");
	l_position = glGetUniformLocation(program, "LightPosition_worldspace");
}
void FEMTest::setGridCellSize()
{
	double length = 0;
	double tn = 0;
	vector<DeformableObject>::iterator object_iter;
	for (object_iter = g_models.begin(); object_iter != g_models.end(); object_iter++)
	{
		length += object_iter->totalLength;
		tn += object_iter->total_tetrahedra;
	}
	g_l = length / (6 * tn);
}
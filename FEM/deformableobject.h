#include <App.h>
#include <vector>
#include <map>
#include "hash.h"
#define EPSILON 0.001f
#define EPS2  EPSILON*EPSILON

struct BoundaryTriangle
{
	GLushort indices[3];
};

struct Tetrahedron {
	int indices[4];			//indices
	float volume;			//volume 
	float plastic[6];		//plasticity values
	glm::vec3 e1, e2, e3;	//edges
	glm::mat3 Re;			//Rotational warp of tetrahedron.
	glm::mat3 Ke[4][4];		//Stiffness element matrix
	glm::vec3 B[4];			//Jacobian of shapefunctions; B=SN =[d/dx  0     0 ][wn 0  0]
	//                                  [0    d/dy   0 ][0 wn  0]
	//									[0     0   d/dz][0  0 wn]
	//									[d/dy d/dx   0 ]
	//									[d/dz  0   d/dx]
	//									[0    d/dz d/dy]
};

struct F_add
{
	int LocalId;
	glm::vec3 f;
};

class DeformableObject
{
public:
	vector<Tetrahedron> tetrahedra;
	vector<BoundaryTriangle> bTriangle;
	vector<glm::vec3> grid;
	vector<float> distanceFild;
	//glm::vec3 dhelp;
	//glm::vec3 dhelp2;
	float gs, gt;
	vector<glm::vec3> gradientFild;
	vector<float> distanceV;
	vector<glm::vec3> gradientV;
	vector<vector<int>> trianglesOfVertices;
	vector<glm::vec3> normaloftriangle;
	vector<glm::vec3> normalofvetex;
	float totalLength = 0;
	float ABmx, ABmy, ABmz, ABMx, ABMy, ABMz;
	float l;//local grid size
	int gi, gj, gk;

	float nu = 0.33f;			//Poisson ratio
	float Y = 500000.0f;		//Young modulus
	float density = 1000.0f;
	float creep = 0.20f;
	float yield = 0.04f;
	float mass_damping = 1.0f;
	float m_max = 0.2f;

	float d15; 
	float d16; 
	float d17; 
	float d18; 

	glm::vec3 D; //Isotropic elasticity matrix D


	glm::vec3 gravity = glm::vec3(0.0f, -9.81f, 0.0f);

	size_t total_points = 0;
	int total_tetrahedra = 0;
	int total_btriangle = 0;
	int selected_index = -1;//////

	vector<glm::vec3> Xi;		//Model coordinates
	vector<glm::vec3> X;		//Current coordinates
	vector<glm::vec3> V;		//Velocity
	vector<float> mass;			//Mass matrix
	vector<glm::vec3> F;		//Force
	vector<bool> IsFixed;		//Fixed point
	glm::mat3 I = glm::mat3(1);	//3x3 identity matrix

	typedef map<int, glm::mat3> matrix_map;
	vector<matrix_map> K_row;
	vector<matrix_map> A_row;
	typedef matrix_map::iterator matrix_iterator;
	vector<glm::vec3> F0;
	vector<glm::vec3> b;

	//For Conjugate Gradient 
	vector<glm::vec3> residual;
	vector<glm::vec3> pre;
	vector<glm::vec3> update;

	float tiny = 1e-010f;       // TODO: Should be user controllable
	float tolerence = 0.001f;
	int i_max = 20;
	bool bUseStiffnessWarping = true;


	GLuint                  meshVBuffer;
	GLuint                  meshNBuffer;
	GLuint                  indiceBuffer;

public:
	DeformableObject();
	~DeformableObject();
	DeformableObject(float x, float y, float z, bool ifFixed);
	DeformableObject(size_t xdim, size_t ydim, size_t zdim, float width, float height, float depth);

	float GetTetraVolume(glm::vec3 e1, glm::vec3 e2, glm::vec3 e3) {
		return  (glm::dot(e1, glm::cross(e2, e3))) / 6.0f;
	}
	void AddBTriangle(int i0, int i1, int i2)
	{
		BoundaryTriangle t;
		t.indices[0] = i0;
		t.indices[1] = i1;
		t.indices[2] = i2;

		bTriangle.push_back(t);
	}
	void AddTetrahedron(int i0, int i1, int i2, int i3) {
		Tetrahedron t;

		t.indices[0] = i0;
		t.indices[1] = i1;
		t.indices[2] = i2;
		t.indices[3] = i3;

		tetrahedra.push_back(t);
	}
	bool IfUseStiffnessWarping(){ return bUseStiffnessWarping; };

	void SetSelectIndex(int i){ selected_index = i; };
	int GetSelectIndex(){ return selected_index; };

	void CalculateK();
	void ReadModelFromFile(float x, float y, float z, bool ifFixed, const char*name);
	void ClearStiffnessAssembly();
	void RecalcMassMatrix();
	void InitializePlastic();
	void OnShutdown();
	void ComputeForces();
	glm::mat3 ortho_normalize(glm::mat3 A);
	void UpdateOrientation();
	void ResetOrientation();
	void StiffnessAssembly();
	void AddPlasticityForce(float dt);
	void DynamicsAssembly(float dt);
	void ConjugateGradientSolver(float dt);
	void UpdatePosition(float dt);

	void GroundCollision();
	void StepPhysics(float dt);

	void RenderModel();

	void                    DistanceFildAndGradientFild();

	//float                   minDistanceBetweenVetexAndTriangle(glm::vec3 v, vector<BoundaryTriangle> *bTriangle, glm::vec3 *dhelp2);
	//float                   DistanceBetweenVT(glm::vec3 v, BoundaryTriangle bt, glm::vec3 *dhelp);
	float                   minDistanceBetweenVetexAndTriangle(glm::vec3 v, vector<BoundaryTriangle> *bTriangle);
	float                   DistanceBetweenVT(glm::vec3 v, BoundaryTriangle bt);

	bool                    IsVinsideObj(glm::vec3 v);

	void                    VetexDistanceAndGradient();
	void                    InterpolateDistanceGradient(int i, glm::vec3 *g);

	void                    GenerateBlocks(size_t xdim, size_t ydim, size_t zdim, float width, float height, float depth);
	void                    Reset();
	void                    firstPass(HashMap *h, int objId);
	void                    secondPass(HashMap *h, int objId);
	void                    F_ADD(int objId);

	void                    CalculateVN();
};
#include "deformableobject.h"
ofstream debug2("debug2.txt");
extern ofstream debug;
extern vector<DeformableObject> g_models;
extern float g_l;//global grid size
extern vector<vector<F_add>> g_F;
DeformableObject::DeformableObject()
{
    d15 = Y / (1.0f + nu) / (1.0f - 2 * nu);
    d16 = (1.0f - nu) * d15;
    d17 = nu * d15;
    d18 = Y / 2 / (1.0f + nu);

	D = glm::vec3(d16, d17, d18); //Isotropic elasticity matrix D

	//ReadModelFromFile("bunny.1");
	ReadModelFromFile(0.0f, 0.0f, 0.0f, false, "bunny_300.1");

	mass.resize(total_points);

	//copy positions to buffer 
	A_row.resize(total_points);
	K_row.resize(total_points);
	b.resize(total_points);
	V.resize(total_points);
	F.resize(total_points);
	F0.resize(total_points);
	residual.resize(total_points);
	update.resize(total_points);
	pre.resize(total_points);

	//fill in V
	memset(&(V[0].x), 0, total_points*sizeof(glm::vec3));


	CalculateK();
	ClearStiffnessAssembly();
	RecalcMassMatrix();
	InitializePlastic();
	//distance file and gradient fild
	for (size_t k = 0; k < tetrahedra.size(); k++)
	{
		glm::vec3 x0 = Xi[tetrahedra[k].indices[0]];
		glm::vec3 x1 = Xi[tetrahedra[k].indices[1]];
		glm::vec3 x2 = Xi[tetrahedra[k].indices[2]];
		glm::vec3 x3 = Xi[tetrahedra[k].indices[3]];

		glm::vec3 e10 = x1 - x0;
		glm::vec3 e20 = x2 - x0;
		glm::vec3 e30 = x3 - x0;
		glm::vec3 e12 = x1 - x2;
		glm::vec3 e23 = x2 - x3;
		glm::vec3 e31 = x3 - x1;
		totalLength += glm::length(e10);
		totalLength += glm::length(e20);
		totalLength += glm::length(e30);
		totalLength += glm::length(e12);
		totalLength += glm::length(e23);
		totalLength += glm::length(e31);

	}
	l = totalLength / (total_tetrahedra * 6);
	DistanceFildAndGradientFild();
	VetexDistanceAndGradient();
	
}
DeformableObject::DeformableObject(float x, float y, float z, bool ifFixed)
{
	d15 = Y / (1.0f + nu) / (1.0f - 2 * nu);
	d16 = (1.0f - nu) * d15;
	d17 = nu * d15;
	d18 = Y / 2 / (1.0f + nu);

	D = glm::vec3(d16, d17, d18); //Isotropic elasticity matrix D

	ReadModelFromFile(x, y, z, ifFixed, "bunny_300.1");

	mass.resize(total_points);

	//copy positions to buffer 
	A_row.resize(total_points);
	K_row.resize(total_points);
	b.resize(total_points);
	V.resize(total_points);
	F.resize(total_points);
	F0.resize(total_points);
	residual.resize(total_points);
	update.resize(total_points);
	pre.resize(total_points);

	//fill in V
	memset(&(V[0].x), 0, total_points*sizeof(glm::vec3));


	CalculateK();
	ClearStiffnessAssembly();
	RecalcMassMatrix();
	InitializePlastic();
	//distance file and gradient fild
	for (size_t k = 0; k < tetrahedra.size(); k++) 
	{
		glm::vec3 x0 = Xi[tetrahedra[k].indices[0]];
		glm::vec3 x1 = Xi[tetrahedra[k].indices[1]];
		glm::vec3 x2 = Xi[tetrahedra[k].indices[2]];
		glm::vec3 x3 = Xi[tetrahedra[k].indices[3]];

		glm::vec3 e10 = x1 - x0;
		glm::vec3 e20 = x2 - x0;
		glm::vec3 e30 = x3 - x0;
		glm::vec3 e12 = x1 - x2;
		glm::vec3 e23 = x2 - x3;
		glm::vec3 e31 = x3 - x1;
		totalLength += glm::length(e10);
		totalLength += glm::length(e20);
		totalLength += glm::length(e30);
		totalLength += glm::length(e12);
		totalLength += glm::length(e23);
		totalLength += glm::length(e31);

	}
	l = totalLength / (total_tetrahedra * 6);
    DistanceFildAndGradientFild();
	VetexDistanceAndGradient();	
}
DeformableObject::DeformableObject(size_t xdim, size_t ydim, size_t zdim, float width, float height, float depth)
{
	d15 = Y / (1.0f + nu) / (1.0f - 2 * nu);
	d16 = (1.0f - nu) * d15;
	d17 = nu * d15;
	d18 = Y / 2 / (1.0f + nu);

	D = glm::vec3(d16, d17, d18); //Isotropic elasticity matrix D

	//ReadModelFromFile("bunny.1");
	GenerateBlocks(xdim, ydim, zdim,  width, height, depth);

	total_tetrahedra = tetrahedra.size();

	total_points = X.size();
	mass.resize(total_points);

	//copy positions to buffer 
	A_row.resize(total_points);
	K_row.resize(total_points);
	b.resize(total_points);
	V.resize(total_points);
	F.resize(total_points);
	F0.resize(total_points);
	residual.resize(total_points);
	update.resize(total_points);
	pre.resize(total_points);

	//fill in V
	memset(&(V[0].x), 0, total_points*sizeof(glm::vec3));

	CalculateK();
	ClearStiffnessAssembly();
	RecalcMassMatrix();
	InitializePlastic();

}
DeformableObject::~DeformableObject()
{
	OnShutdown();
}

void DeformableObject::CalculateK() {

	for (size_t k = 0; k<tetrahedra.size(); k++) {

		glm::vec3 x0 = Xi[tetrahedra[k].indices[0]];
		glm::vec3 x1 = Xi[tetrahedra[k].indices[1]];
		glm::vec3 x2 = Xi[tetrahedra[k].indices[2]];
		glm::vec3 x3 = Xi[tetrahedra[k].indices[3]];

		//For this check page no.: 344-346 of Kenny Erleben's book Physics based Animation
		//Eq. 10.30(a-c)
		glm::vec3 e10 = x1 - x0;
		glm::vec3 e20 = x2 - x0;
		glm::vec3 e30 = x3 - x0;

		tetrahedra[k].e1 = e10;
		tetrahedra[k].e2 = e20;
		tetrahedra[k].e3 = e30;

		tetrahedra[k].volume = GetTetraVolume(e10, e20, e30);

		//Eq. 10.32
		glm::mat3 E = glm::mat3(e10.x, e10.y, e10.z,
			e20.x, e20.y, e20.z,
			e30.x, e30.y, e30.z);
		float detE = glm::determinant(E);
		float invDetE = 1.0f / detE;

		//Eq. 10.40 (a) & Eq. 10.42 (a)
		//Shape function derivatives wrt x,y,z
		// d/dx N0
		float invE10 = (e20.z*e30.y - e20.y*e30.z)*invDetE;
		float invE20 = (e10.y*e30.z - e10.z*e30.y)*invDetE;
		float invE30 = (e10.z*e20.y - e10.y*e20.z)*invDetE;
		float invE00 = -invE10 - invE20 - invE30;

		//Eq. 10.40 (b) & Eq. 10.42 (b)
		// d/dy N0
		float invE11 = (e20.x*e30.z - e20.z*e30.x)*invDetE;
		float invE21 = (e10.z*e30.x - e10.x*e30.z)*invDetE;
		float invE31 = (e10.x*e20.z - e10.z*e20.x)*invDetE;
		float invE01 = -invE11 - invE21 - invE31;

		//Eq. 10.40 (c) & Eq. 10.42 (c)
		// d/dz N0
		float invE12 = (e20.y*e30.x - e20.x*e30.y)*invDetE;
		float invE22 = (e10.x*e30.y - e10.y*e30.x)*invDetE;
		float invE32 = (e10.y*e20.x - e10.x*e20.y)*invDetE;
		float invE02 = -invE12 - invE22 - invE32;

		glm::vec3 B[4];
		//Eq. 10.43 
		//Bn ~ [bn cn dn]^T
		// bn = d/dx N0 = [ invE00 invE10 invE20 invE30 ]
		// cn = d/dy N0 = [ invE01 invE11 invE21 invE31 ]
		// dn = d/dz N0 = [ invE02 invE12 invE22 invE32 ]
		tetrahedra[k].B[0] = glm::vec3(invE00, invE01, invE02);
		tetrahedra[k].B[1] = glm::vec3(invE10, invE11, invE12);
		tetrahedra[k].B[2] = glm::vec3(invE20, invE21, invE22);
		tetrahedra[k].B[3] = glm::vec3(invE30, invE31, invE32);

		for (int i = 0; i<4; i++) {
			for (int j = 0; j<4; j++) {
				glm::mat3& Ke = tetrahedra[k].Ke[i][j];
				float d19 = tetrahedra[k].B[i].x;
				float d20 = tetrahedra[k].B[i].y;
				float d21 = tetrahedra[k].B[i].z;
				float d22 = tetrahedra[k].B[j].x;
				float d23 = tetrahedra[k].B[j].y;
				float d24 = tetrahedra[k].B[j].z;
				Ke[0][0] = d16 * d19 * d22 + d18 * (d20 * d23 + d21 * d24);
				Ke[0][1] = d17 * d19 * d23 + d18 * (d20 * d22);
				Ke[0][2] = d17 * d19 * d24 + d18 * (d21 * d22);

				Ke[1][0] = d17 * d20 * d22 + d18 * (d19 * d23);
				Ke[1][1] = d16 * d20 * d23 + d18 * (d19 * d22 + d21 * d24);
				Ke[1][2] = d17 * d20 * d24 + d18 * (d21 * d23);

				Ke[2][0] = d17 * d21 * d22 + d18 * (d19 * d24);
				Ke[2][1] = d17 * d21 * d23 + d18 * (d20 * d24);
				Ke[2][2] = d16 * d21 * d24 + d18 * (d20 * d23 + d19 * d22);

				Ke *= tetrahedra[k].volume;
			}
		}
	}
}
void DeformableObject::ReadModelFromFile(float x, float y, float z, bool ifFixed, const char *filename)
{
	string fileName = filename;
	string suffix = ".node";
	string suffix2 = ".ele";
	string suffix3 = ".face";
	fileName += suffix;
	ifstream nodeFile(fileName.c_str(), ios::in);
	if (!nodeFile.is_open())
	{
		cout << fileName << " open fail." << endl;
		return;
	}

	nodeFile >> total_points;
	X.resize(total_points);
	Xi.resize(total_points);
	IsFixed.resize(total_points);
	trianglesOfVertices.resize(total_points);
	normalofvetex.resize(total_points);
	normaloftriangle.resize(total_points);
	int noUse;
	nodeFile >> noUse >> noUse >> noUse;

	glm::vec3 nodePosition;
	float miny = 10;
	for (int i = 0; i < total_points; ++i)
	{
		nodeFile >> noUse;
		nodeFile >> nodePosition.x >> nodePosition.y >> nodePosition.z;
		//debug << nodePosition.x << " " << nodePosition.y << " " << nodePosition.z << endl;
		nodePosition.x += x;
		nodePosition.y += y;
		nodePosition.z += z;
		X[i] = nodePosition;
		Xi[i] = X[i];
		//Make the first few points fixed

		if (ifFixed)
		if (Xi[i].y < 0.04)
		{
			miny = glm::min(miny, Xi[i].y);
			IsFixed[i] = true;
		}
		else
		{
			IsFixed[i] = false;
		}
	}
	//move to floor
	//float k = X[ind].y;

	if (ifFixed)
	for (size_t i = 0; i<total_points; i++) {
		X[i].y -= miny;
		Xi[i].y = X[i].y;
	}
	nodeFile.close();

	fileName = filename;
	fileName += suffix2;

	ifstream eleFile(fileName.c_str(), ios::in);
	if (!eleFile.is_open())
	{
		cout << fileName << " open fail." << endl;
		return;
	}

	eleFile >> total_tetrahedra;
	eleFile >> noUse >> noUse;

	for (int i = 0; i < total_tetrahedra; ++i)
	{
		int p0, p1, p2, p3;
		eleFile >> noUse;
		eleFile >> p0 >> p1 >> p2 >> p3;
		//debug << p0 << " " << p1 << " " << p2 << " " << p3 <<endl;
		AddTetrahedron(p0, p1, p2, p3);
	}
	eleFile.close();


	fileName = filename;
	fileName += suffix3;

	distanceV.resize(total_points);
	gradientV.resize(total_points);
	for (int i = 0; i < total_points; ++i)
	{
		distanceV[i] = 1;
	}

	ifstream faceFile(fileName.c_str(), ios::in);
	if (!faceFile.is_open())
	{
	cout << fileName << " open fail." << endl;
	return;
	}
	faceFile >> total_btriangle;
	faceFile >> noUse;
	normaloftriangle.resize(total_btriangle);
	for (int i = 0; i < total_btriangle; ++i)
	{
		int p0, p1, p2;
		faceFile >> noUse;
		faceFile >> p0 >> p1 >> p2;
		distanceV[p0] = 0;
		distanceV[p1] = 0;
		distanceV[p2] = 0;
		trianglesOfVertices[p0].push_back(i);
		trianglesOfVertices[p1].push_back(i);
		trianglesOfVertices[p2].push_back(i);
		AddBTriangle(p0, p1, p2);
	}
	faceFile.close();
}


void DeformableObject::RecalcMassMatrix() {
	//This is a lumped mass matrix
	//Based on Eq. 10.106 and pseudocode in Fig. 10.9 on page 358
	for (size_t i = 0; i<total_points; i++) {
		if (IsFixed[i])
			mass[i] = std::numeric_limits<float>::max();
		else
			mass[i] = 1.0f / total_points;
	}

	for (int i = 0; i<total_tetrahedra; i++) {
		float m = (density*tetrahedra[i].volume)* 0.25f;
		mass[tetrahedra[i].indices[0]] += m;
		mass[tetrahedra[i].indices[1]] += m;
		mass[tetrahedra[i].indices[2]] += m;
		mass[tetrahedra[i].indices[3]] += m;
	}
}

void DeformableObject::InitializePlastic() {
	for (size_t i = 0; i<tetrahedra.size(); i++) {
		for (int j = 0; j<6; j++)
			tetrahedra[i].plastic[j] = 0;
	}
}

void DeformableObject::ClearStiffnessAssembly() {
	for (size_t k = 0; k<total_points; k++) {
		F0[k].x = 0.0f;
		F0[k].y = 0.0f;
		F0[k].z = 0.0f;

		for (matrix_iterator Kij = K_row[k].begin(); Kij != K_row[k].end(); ++Kij)
			Kij->second = glm::mat3(0);
	}
}


void DeformableObject::OnShutdown() {
	X.clear();
	Xi.clear();
	V.clear();
	mass.clear();
	F.clear();
	IsFixed.clear();
	tetrahedra.clear();
	K_row.clear();
	A_row.clear();
	F0.clear();
	b.clear();
	residual.clear();
	pre.clear();
	update.clear();
}



void DeformableObject::ComputeForces() {
	size_t i = 0;
	for (i = 0; i<total_points; i++) {

		//add gravity force only for non-fixed points
		F[i] += gravity*mass[i];
	}
}


glm::mat3 DeformableObject::ortho_normalize(glm::mat3 A) {
	glm::vec3 row0(A[0][0], A[0][1], A[0][2]);
	glm::vec3 row1(A[1][0], A[1][1], A[1][2]);
	glm::vec3 row2(A[2][0], A[2][1], A[2][2]);

	float L0 = glm::length(row0);
	if (L0)
		row0 /= L0;

	row1 -= row0 * glm::dot(row0, row1);
	float L1 = glm::length(row1);
	if (L1)
		row1 /= L1;

	row2 = glm::cross(row0, row1);

	return glm::mat3(row0,
		row1,
		row2);
}

void DeformableObject::UpdateOrientation() {
	for (int k = 0; k<total_tetrahedra; k++) {
		//Based on description on page 362-364 
		float div6V = 1.0f / tetrahedra[k].volume*6.0f;

		int i0 = tetrahedra[k].indices[0];
		int i1 = tetrahedra[k].indices[1];
		int i2 = tetrahedra[k].indices[2];
		int i3 = tetrahedra[k].indices[3];

		glm::vec3 p0 = X[i0];
		glm::vec3 p1 = X[i1];
		glm::vec3 p2 = X[i2];
		glm::vec3 p3 = X[i3];

		glm::vec3 e1 = p1 - p0;
		glm::vec3 e2 = p2 - p0;
		glm::vec3 e3 = p3 - p0;

		//Eq. 10.129 on page 363 & Eq. 10.131 page 364
		//n1,n2,n3 approximate Einv
		glm::vec3 n1 = glm::cross(e2, e3) * div6V;
		glm::vec3 n2 = glm::cross(e3, e1) * div6V;
		glm::vec3 n3 = glm::cross(e1, e2) * div6V;

		//Now get the rotation matrix from the initial undeformed (model/material coordinates)
		//We get the precomputed edge values
		e1 = tetrahedra[k].e1;
		e2 = tetrahedra[k].e2;
		e3 = tetrahedra[k].e3;

		//Based on Eq. 10.133		

		tetrahedra[k].Re[0][0] = e1.x * n1.x + e2.x * n2.x + e3.x * n3.x;
		tetrahedra[k].Re[0][1] = e1.x * n1.y + e2.x * n2.y + e3.x * n3.y;
		tetrahedra[k].Re[0][2] = e1.x * n1.z + e2.x * n2.z + e3.x * n3.z;

		tetrahedra[k].Re[1][0] = e1.y * n1.x + e2.y * n2.x + e3.y * n3.x;
		tetrahedra[k].Re[1][1] = e1.y * n1.y + e2.y * n2.y + e3.y * n3.y;
		tetrahedra[k].Re[1][2] = e1.y * n1.z + e2.y * n2.z + e3.y * n3.z;

		tetrahedra[k].Re[2][0] = e1.z * n1.x + e2.z * n2.x + e3.z * n3.x;
		tetrahedra[k].Re[2][1] = e1.z * n1.y + e2.z * n2.y + e3.z * n3.y;
		tetrahedra[k].Re[2][2] = e1.z * n1.z + e2.z * n2.z + e3.z * n3.z;

		tetrahedra[k].Re = ortho_normalize(tetrahedra[k].Re);

	}
}
void DeformableObject::ResetOrientation() {
	for (int k = 0; k<total_tetrahedra; k++) {
		tetrahedra[k].Re = I;
	}
}

void DeformableObject::StiffnessAssembly() {

	for (int k = 0; k<total_tetrahedra; k++) {
		glm::mat3 Re = tetrahedra[k].Re;
		glm::mat3 ReT = glm::transpose(Re);


		for (int i = 0; i < 4; ++i) {
			//Based on pseudocode given in Fig. 10.11 on page 361
			glm::vec3 f = glm::vec3(0.0f, 0.0f, 0.0f);
			for (int j = 0; j < 4; ++j) {
				glm::mat3 tmpKe = tetrahedra[k].Ke[i][j];
				glm::vec3 x0 = Xi[tetrahedra[k].indices[j]];
				glm::vec3 prod = glm::vec3(tmpKe[0][0] * x0.x + tmpKe[0][1] * x0.y + tmpKe[0][2] * x0.z, //tmpKe*x0;
					tmpKe[1][0] * x0.x + tmpKe[1][1] * x0.y + tmpKe[1][2] * x0.z,
					tmpKe[2][0] * x0.x + tmpKe[2][1] * x0.y + tmpKe[2][2] * x0.z);
				f += prod;
				if (j >= i) {
					//Based on pseudocode given in Fig. 10.12 on page 361
					glm::mat3 tmp = Re*tmpKe*ReT;
					int index = tetrahedra[k].indices[i];

					K_row[index][tetrahedra[k].indices[j]] += (tmp);

					if (j > i) {
						index = tetrahedra[k].indices[j];
						K_row[index][tetrahedra[k].indices[i]] += (glm::transpose(tmp));
					}
				}

			}
			int idx = tetrahedra[k].indices[i];
			F0[idx] -= Re*f;
		}
	}
}

void DeformableObject::AddPlasticityForce(float dt) {
	for (int k = 0; k<total_tetrahedra; k++) {
		float e_total[6];
		float e_elastic[6];
		for (int i = 0; i<6; ++i)
			e_elastic[i] = e_total[i] = 0;

		//--- Compute total strain: e_total  = Be (Re^{-1} x - x0)
		for (unsigned int j = 0; j<4; ++j) {

			glm::vec3 x_j = X[tetrahedra[k].indices[j]];
			glm::vec3 x0_j = Xi[tetrahedra[k].indices[j]];
			glm::mat3 ReT = glm::transpose(tetrahedra[k].Re);
			glm::vec3 prod = glm::vec3(ReT[0][0] * x_j.x + ReT[0][1] * x_j.y + ReT[0][2] * x_j.z, //tmpKe*x0;
				ReT[1][0] * x_j.x + ReT[1][1] * x_j.y + ReT[1][2] * x_j.z,
				ReT[2][0] * x_j.x + ReT[2][1] * x_j.y + ReT[2][2] * x_j.z);

			glm::vec3 tmp = prod - x0_j;

			//B contains Jacobian of shape funcs. B=SN
			float bj = tetrahedra[k].B[j].x;
			float cj = tetrahedra[k].B[j].y;
			float dj = tetrahedra[k].B[j].z;

			e_total[0] += bj*tmp.x;
			e_total[1] += cj*tmp.y;
			e_total[2] += dj*tmp.z;
			e_total[3] += cj*tmp.x + bj*tmp.y;
			e_total[4] += dj*tmp.x + bj*tmp.z;
			e_total[5] += dj*tmp.y + cj*tmp.z;
		}

		//--- Compute elastic strain
		for (int i = 0; i<6; ++i)
			e_elastic[i] = e_total[i] - tetrahedra[k].plastic[i];

		//--- if elastic strain exceeds c_yield then it is added to plastic strain by c_creep
		float norm_elastic = 0;
		for (int i = 0; i<6; ++i)
			norm_elastic += e_elastic[i] * e_elastic[i];
		norm_elastic = sqrt(norm_elastic);
		if (norm_elastic > yield) {
			float amount = dt*glm::min(creep, (1.0f / dt));  //--- make sure creep do not exceed 1/dt
			for (int i = 0; i<6; ++i)
				tetrahedra[k].plastic[i] += amount*e_elastic[i];
		}

		//--- if plastic strain exceeds c_max then it is clamped to maximum magnitude
		float norm_plastic = 0;
		for (int i = 0; i<6; ++i)
			norm_plastic += tetrahedra[k].plastic[i] * tetrahedra[k].plastic[i];
		norm_plastic = sqrt(norm_plastic);

		if (norm_plastic > m_max) {
			float scale = m_max / norm_plastic;
			for (int i = 0; i<6; ++i)
				tetrahedra[k].plastic[i] *= scale;
		}

		for (size_t n = 0; n<4; ++n) {
			float* e_plastic = tetrahedra[k].plastic;
			//bn, cn and dn are the shape function derivative wrt. x,y and z axis
			//These were calculated in CalculateK function

			//Eq. 10.140(a) & (b) on page 365
			float bn = tetrahedra[k].B[n].x;
			float cn = tetrahedra[k].B[n].y;
			float dn = tetrahedra[k].B[n].z;
			float D0 = D.x;
			float D1 = D.y;
			float D2 = D.z;
			glm::vec3 f = glm::vec3(0);

			float  bnD0 = bn*D0;
			float  bnD1 = bn*D1;
			float  bnD2 = bn*D2;
			float  cnD0 = cn*D0;
			float  cnD1 = cn*D1;
			float  cnD2 = cn*D2;
			float  dnD0 = dn*D0;
			float  dnD1 = dn*D1;
			float  dnD2 = dn*D2;

			//Eq. 10.141 on page 365
			f.x = bnD0*e_plastic[0] + bnD1*e_plastic[1] + bnD1*e_plastic[2] + cnD2*e_plastic[3] + dnD2*e_plastic[4];
			f.y = cnD1*e_plastic[0] + cnD0*e_plastic[1] + cnD1*e_plastic[2] + bnD2*e_plastic[3] + +dnD2*e_plastic[5];
			f.z = dnD1*e_plastic[0] + dnD1*e_plastic[1] + dnD0*e_plastic[2] + bnD2*e_plastic[4] + cnD2*e_plastic[5];

			f *= tetrahedra[k].volume;
			int idx = tetrahedra[k].indices[n];
			F[idx] += tetrahedra[k].Re*f;
		}
	}
}
void DeformableObject::DynamicsAssembly(float dt) {
	float dt2 = dt*dt;

	for (size_t k = 0; k<total_points; k++) {

		float m_i = mass[k];
		b[k] = glm::vec3(0.0, 0.0, 0.0);
		//printf("Idx: %3d\n",k);

		matrix_map tmp = K_row[k];
		matrix_iterator Kbegin = tmp.begin();
		matrix_iterator Kend = tmp.end();
		for (matrix_iterator K = Kbegin; K != Kend; ++K)
		{
			unsigned int j = K->first;
			glm::mat3 K_ij = K->second;
			glm::vec3 x_j = X[j];
			glm::mat3& A_ij = A_row[k][j];

			A_ij = K_ij * (dt2);
			glm::vec3 prod = glm::vec3(K_ij[0][0] * x_j.x + K_ij[0][1] * x_j.y + K_ij[0][2] * x_j.z,
				K_ij[1][0] * x_j.x + K_ij[1][1] * x_j.y + K_ij[1][2] * x_j.z,
				K_ij[2][0] * x_j.x + K_ij[2][1] * x_j.y + K_ij[2][2] * x_j.z);

			b[k] -= prod;//K_ij * x_j;

			if (k == j)
			{
				float c_i = mass_damping*m_i;
				float tmp = m_i + dt*c_i;
				A_ij[0][0] += tmp;
				A_ij[1][1] += tmp;
				A_ij[2][2] += tmp;
			}
		}

		b[k] -= F0[k];
		b[k] += F[k];
		b[k] *= dt;
		b[k] += V[k] * m_i;
	}
}
void DeformableObject::ConjugateGradientSolver(float dt) {


	for (size_t k = 0; k<total_points; k++) {
		if (IsFixed[k])
			continue;
		residual[k] = b[k];

		matrix_iterator Abegin = A_row[k].begin();
		matrix_iterator Aend = A_row[k].end();
		for (matrix_iterator A = Abegin; A != Aend; ++A)
		{
			unsigned int j = A->first;
			glm::mat3& A_ij = A->second;
			float v_jx = V[j].x;
			float v_jy = V[j].y;
			float v_jz = V[j].z;
			glm::vec3 prod = glm::vec3(A_ij[0][0] * v_jx + A_ij[0][1] * v_jy + A_ij[0][2] * v_jz, //A_ij * pre[j]
				A_ij[1][0] * v_jx + A_ij[1][1] * v_jy + A_ij[1][2] * v_jz,
				A_ij[2][0] * v_jx + A_ij[2][1] * v_jy + A_ij[2][2] * v_jz);
			residual[k] -= prod;//  A_ij * v_j;

		}
		pre[k] = residual[k];
	}

	for (int i = 0; i<i_max; i++) {
		float d = 0;
		float d2 = 0;

		for (size_t k = 0; k<total_points; k++) {

			if (IsFixed[k])
				continue;

			update[k].x = 0; update[k].y = 0; update[k].z = 0;

			matrix_iterator Abegin = A_row[k].begin();
			matrix_iterator Aend = A_row[k].end();
			for (matrix_iterator A = Abegin; A != Aend; ++A) {
				unsigned int j = A->first;
				glm::mat3& A_ij = A->second;
				float prevx = pre[j].x;
				float prevy = pre[j].y;
				float prevz = pre[j].z;
				glm::vec3 prod = glm::vec3(A_ij[0][0] * prevx + A_ij[0][1] * prevy + A_ij[0][2] * prevz, //A_ij * pre[j]
					A_ij[1][0] * prevx + A_ij[1][1] * prevy + A_ij[1][2] * prevz,
					A_ij[2][0] * prevx + A_ij[2][1] * prevy + A_ij[2][2] * prevz);
				update[k] += prod;//A_ij*pre[j];

			}
			d += glm::dot(residual[k], residual[k]);
			d2 += glm::dot(pre[k], update[k]);
		}

		if (fabs(d2)<tiny)
			d2 = tiny;

		float d3 = d / d2;
		float d1 = 0;


		for (size_t k = 0; k<total_points; k++) {
			if (IsFixed[k])
				continue;

			V[k] += pre[k] * d3;
			residual[k] -= update[k] * d3;
			d1 += glm::dot(residual[k], residual[k]);
		}

		if (i >= i_max && d1 < tolerence)
			break;

		if (fabs(d)<tiny)
			d = tiny;

		float d4 = d1 / d;

		for (size_t k = 0; k<total_points; k++) {
			if (IsFixed[k])
				continue;
			pre[k] = residual[k] + pre[k] * d4;
		}
	}
}
void DeformableObject::UpdatePosition(float dt) {
	for (size_t k = 0; k<total_points; k++) {
		if (IsFixed[k])
			continue;
		X[k] += float(dt)*V[k];
	}
}


void DeformableObject::GroundCollision()
{
	for (size_t i = 0; i<total_points; i++) {
		if (X[i].y < 0) //collision with ground
		{
			glm::vec3 f;
			f = glm::vec3(0.0f, -3000.0f, 0.0f) * X[i].y;
			F[i] += f;
			//X[i].y = 0;
		}
	}
}

void DeformableObject::StepPhysics(float dt) {
	ComputeForces();

	GroundCollision();

	ClearStiffnessAssembly();

	if (bUseStiffnessWarping)
		UpdateOrientation();
	else
		ResetOrientation();

	StiffnessAssembly();

	AddPlasticityForce(dt);

	DynamicsAssembly(dt);

	ConjugateGradientSolver(dt);
	//EigenSolve();

	UpdatePosition(dt);
	//GroundCollision();
}

void DeformableObject::RenderModel()
{
	/*
	glColor3f(0.75, 0.75, 0.75);
	glBegin(GL_LINES);

	for (int i = 0; i<total_tetrahedra; i++) {
		int i0 = tetrahedra[i].indices[0];
		int i1 = tetrahedra[i].indices[1];
		int i2 = tetrahedra[i].indices[2];
		int i3 = tetrahedra[i].indices[3];
		glm::vec3 p1 = X[i0];
		glm::vec3 p2 = X[i1];
		glm::vec3 p3 = X[i2];
		glm::vec3 p4 = X[i3];

		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p3.x, p3.y, p3.z);

		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p3.x, p3.y, p3.z);

		glVertex3f(p2.x, p2.y, p2.z);		glVertex3f(p3.x, p3.y, p3.z);
	}
	glEnd();
*/

	//glColor3f(0.75, 0.75, 0.75);
	glGenBuffers(1, &meshVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		X.size() * sizeof(glm::vec3),
		&(X[0]),
		GL_STATIC_DRAW);
	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

	glGenBuffers(1, &meshNBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, meshNBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		normalofvetex.size() * sizeof(glm::vec3),
		&(normalofvetex[0]),
		GL_STATIC_DRAW);
	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, meshNBuffer);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);


	glGenBuffers(1, &indiceBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indiceBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, bTriangle.size() * sizeof(BoundaryTriangle), &(bTriangle[0]), GL_STATIC_DRAW);
	
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indiceBuffer);

	// Draw the triangles !
	glDrawElements(
		GL_TRIANGLES,      // mode
		bTriangle.size() * 3,    // count
		GL_UNSIGNED_SHORT,   // type
		(void*)0           // element array buffer offset
		);

	/*
	glBegin(GL_TRIANGLES);

	for (int i = 0; i<bTriangle.size(); i++) {
		int i0 = bTriangle[i].indices[0];
		int i1 = bTriangle[i].indices[1];
		int i2 = bTriangle[i].indices[2];

		glm::vec3 p1 = X[i0];
		glm::vec3 p2 = X[i1];
		glm::vec3 p3 = X[i2];


		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p3.x, p3.y, p3.z);
	}
	glEnd();*/
	/*
	//draw points	
	glBegin(GL_POINTS);
	for (int i = 0; i<total_points; i++) {
		glm::vec3 p = X[i];
		int is = (i == selected_index);
		glColor3f((float)!is, (float)is, (float)is);
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();*/
}

void DeformableObject::DistanceFildAndGradientFild()
{
	float mx, my, mz, Mx, My, Mz;
	mx = X[0].x;   Mx = X[0].x;
	my = X[0].y;   My = X[0].y;
	mz = X[0].z;   Mz = X[0].z;
	for (int i = 0; i < total_points; ++i)
	{
		mx = glm::min(X[i].x, mx);
		my = glm::min(X[i].y, my);
		mz = glm::min(X[i].z, mz);
		Mx = glm::max(X[i].x, Mx);
		My = glm::max(X[i].y, My);
		Mz = glm::max(X[i].z, Mz);
	}
	//debug << mx << " " << my << " " << mz << " " << Mx << " " << My << " " << Mz << endl;


	mx -= l;
	my -= l;
	mz -= l;
	Mx += l;
	My += l;
	Mz += l;
	Mx += l;
	My += l;
	Mz += l;
	ABmx = mx;
	ABmy = my;
	ABmz = mz;
	ABMx = Mx;
	ABMy = My;
	ABMz = Mz;
	//debug << mx << " " << my << " " << mz << " " << Mx << " " << My << " " << Mz << endl;
	float li = Mx - mx;
	float lj = My - my;
	float lk = Mz - mz;
	int xi = li / l + 1;
	int yj = lj / l + 1;
	int zk = lk / l + 1;

	gi = xi;
	gj = yj;
	gk = zk;
	//debug << l << " " << gi << " " << gj << " " << gk << endl;
	//cout << li <<"  "<<l<<"  "<< xi << " " << yj << " " << zk << endl;
	grid.resize(xi*yj*zk);
	distanceFild.resize(xi*yj*zk);
	gradientFild.resize(xi*yj*zk);

	for (int i = 0; i < xi; ++i)
		for (int j = 0; j < yj; ++j)
			for (int k = 0; k < zk; ++k)
			{
		int index = i * yj * zk + j * zk + k;
		grid[index] = glm::vec3(mx + (i * l), my + (j * l), mz + (k * l));
		distanceFild[index] = - minDistanceBetweenVetexAndTriangle(grid[index], &bTriangle);
		if (IsVinsideObj(grid[index]))
		{ 
			distanceFild[index] = -distanceFild[index];
		}

		//distanceFildV[index] = dhelp2;
		//debug << distanceFild[index] << endl;
			}
	//debug << "distanceFild over." << endl;

	for (int i = 1; i < xi - 1; ++i)
		for (int j = 1; j < yj - 1; ++j)
			for (int k = 1; k < zk - 1; ++k)
			{
		int di = i * yj * zk;
		int dj = j * zk;
		int dk = k;
		int dip1 = di + yj * zk;
		int di_1 = di - yj * zk;
		int djp1 = dj + zk;
		int dj_1 = dj - zk;
		int dkp1 = dk + 1;
		int dk_1 = dk - 1;
		gradientFild[di + dj + dk] = -glm::vec3(
			distanceFild[dip1 + dj + dk] - distanceFild[di_1 + dj + dk], 
			distanceFild[di + djp1 + dk] - distanceFild[di_1 + dj_1 + dk], 
			distanceFild[di + dj + dkp1] - distanceFild[di + dj + dk_1]);
			}
}
bool  DeformableObject::IsVinsideObj(glm::vec3 v)
{
	for (int i = 0; i < total_tetrahedra; ++i)
	{
		int n1 = tetrahedra[i].indices[0];
		int n2 = tetrahedra[i].indices[1];
		int n3 = tetrahedra[i].indices[2];
		int n4 = tetrahedra[i].indices[3];
		glm::vec3 v0 = X[n1];
		glm::vec3 v1 = X[n2];
		glm::vec3 v2 = X[n3];
		glm::vec3 v3 = X[n4];

		glm::mat3 A;
		glm::vec3 e1, e2, e3;
		e1 = v1 - v0;
		e2 = v2 - v0;
		e3 = v3 - v0;

		A[0][0] = e1.x;    A[0][1] = e2.x;    A[0][2] = e3.x;
		A[1][0] = e1.y;    A[1][1] = e2.y;    A[1][2] = e3.y;
		A[2][0] = e1.z;    A[2][1] = e2.z;    A[2][2] = e3.z;

		A = glm::inverse(A);
		glm::vec3 beta = v - v0;
		glm::vec3 b;

		b.x = A[0][0] * beta.x + A[0][1] * beta.y + A[0][2] * beta.z;
		b.y = A[1][0] * beta.x + A[1][1] * beta.y + A[1][2] * beta.z;
		b.z = A[2][0] * beta.x + A[2][1] * beta.y + A[2][2] * beta.z;


		if (b.x >= 0 && b.y >= 0 && b.z >= 0 && (b.x + b.y + b.z) <= 1)
			return true;
	}
	return false;
}
float DeformableObject::minDistanceBetweenVetexAndTriangle(glm::vec3 v, vector<BoundaryTriangle> *bTriangle)
//float DeformableObject::minDistanceBetweenVetexAndTriangle(glm::vec3 v, vector<BoundaryTriangle> *bTriangle, glm::vec3 *dhelp2)
{
	float distance = 100000;
	for (int i = 0; i < total_btriangle; ++i)
	{
		float d = DistanceBetweenVT(v, (*bTriangle)[i]);
		if (d < distance)
		{
			distance = d;
			//*dhelp2 = dhelp;
		}
	}
	return distance;
}
float DeformableObject::DistanceBetweenVT(glm::vec3 v, BoundaryTriangle bt)
//float DeformableObject::DistanceBetweenVT(glm::vec3 v, BoundaryTriangle bt, glm::vec3 *dhelp)
{//geometric tools for computer graphics p.275
	glm::vec3 v0, v1, v2, dv;
	float a, b, c, d, e, f;
	v0 = X[bt.indices[0]];
	v1 = X[bt.indices[1]];
	v2 = X[bt.indices[2]];
	v1 = v1 - v0;//e0
	v2 = v2 - v0;//e1
	dv = v0 - v;

	a = glm::dot(v1, v1);
	b = glm::dot(v1, v2);
	c = glm::dot(v2, v2);
	d = glm::dot(v1, dv);
	e = glm::dot(v2, dv);
	f = glm::dot(dv, dv);

	float det, s, t;
	det = a * c - b * b;
	s = b * e - c * d;
	t = b * d - a * e;

	if (s >= 0 && s <= det && t >= 0 && t <= det && (s + t) <= det)//p` in the triangle.
	{
		float invDet = 1 / det;
		s *= invDet;
		t *= invDet;
		//*dhelp = v0 + s*v1 + t*v2;
		return sqrt(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f);
	}
	else
		if ((s + t) <= det)
	    {
		if (s < 0)
		{
			if (t < 0)
			{//region 4
				if (d > 0)
				{//on t = 0;
					t = 0;
					s = (d >= 0 ? 0 : (-d >= a ? 1 : -d / a));
					//*dhelp = v0 + s*v1 + t*v2;
					return sqrt(a*s*s + 2 * d*s + f);
				}
				else
				{
					s = 0;
					t = (e >= 0 ? 0 : (-e >= c ? 1 : -e / c));
					//*dhelp = v0 + s*v1 + t*v2;
					return sqrt(c*t*t + 2 * e*t + f);
				}
			}
			else
			{//region 3
				s = 0;
				t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
				//*dhelp = v0 + s*v1 + t*v2;
				return sqrt(c*t*t + 2 * e*t + f);
			}
		}
		else
		{//region 5
			t = 0;
			s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
			//*dhelp = v0 + s*v1 + t*v2;
			return sqrt(a*s*s + 2 * d*s + 2 * e*t + f);
		}
	    }
		else
		{
			if (s < 0)
			{//region 2
				float tmp0 = b + d;
				float tmp1 = c + e;
				if (tmp1 > tmp0)
				{//min on edge s + t = 1
					float numer = tmp1 - tmp0;
					float denom = a - 2 * b + c;
					s = (numer >= denom ? 1 : numer/denom);
					t = 1 - s;
					//*dhelp = v0 + s*v1 + t*v2;
					return sqrt(a*s*s + 2 * b*s*t + c*t*t + 2 * d*s + 2 * e*t + f);
				}
				else
				{
					s = 0;
					t = (tmp1 < 0 ? 1 : (e >= 0 ? 0 : -e/c));
					//*dhelp = v0 + s*v1 + t*v2;
					return sqrt(c*t*t + 2 * e*t + f);
				}
			}
			else if (t < 0)
			{//region 6
				float tmp0 = b + e;
				float tmp1 = a + d;
				if (tmp1 > tmp0)
				{//min on edge s + t = 1
					float numer = c + e - b - d;
					float denom = a - 2 * b + c;
					s = (numer >= denom ? 1 : numer / denom);
					t = 1 - s;
					//*dhelp = v0 + s*v1 + t*v2;
					return sqrt(a*s*s + 2 * b*s*t + c*t*t + 2 * d*s + 2 * e*t + f);
				}
				else
				{// or on 
					t = 0;
					s = (tmp1 < 0 ? 1 : (d >= 0 ? 0 : -d / a));
					//*dhelp = v0 + s*v1 + t*v2;
					return sqrt(a*s*s + 2 * d*s + f);
				}
			}
			else
			{//region 1
				float numer = c + d - b - d;
				if (numer <= 0)
				{
					s = 0;
				}
				else
				{
					float denom = a - 2 * b + c;
					s = (numer >= denom ? 1 : numer / denom);
				}
				t = 1 - s;
				//*dhelp = v0 + s*v1 + t*v2;
				return sqrt(a*s*s + 2 * b*s*t + c*t*t + 2 * d*s + 2 * e*t + f);
			}
		}
}
void DeformableObject::VetexDistanceAndGradient()
{
	for (int i = 0; i < total_points; ++i)
	{
		if (distanceV[i])
			distanceV[i] = minDistanceBetweenVetexAndTriangle(X[i], &bTriangle);
		InterpolateDistanceGradient(i, &gradientV[i]);
	}
}
void DeformableObject::InterpolateDistanceGradient(int ind, glm::vec3 *g)
{
	float li, lj, lk;
	li = X[ind].x - ABmx;
	lj = X[ind].y - ABmy;
	lk = X[ind].z - ABmz;
	int i, j, k;
	i = li / l;
	j = lj / l;
	k = lk / l;

	glm::vec3 a;
	glm::vec3 ga, gb, gc, gd, ge, gf, gg, gh;
	float da, db, dc, dd, de, df, dg, dh;
	a = glm::vec3(ABmx + i * l, ABmy + j * l, ABmz + k * l);
	//b = glm::vec3(ABmx + i * l + l, ABmy + j * l, ABmz + k * l);
	//c = glm::vec3(ABmx + i * l + l, ABmy + j * l, ABmz + k * l + l);
	//d = glm::vec3(ABmx + i * l, ABmy + j * l, ABmz + k * l + l);

	//e = glm::vec3(ABmx + i * l, ABmy + j * l + l, ABmz + k * l);
	//f = glm::vec3(ABmx + i * l + l, ABmy + j * l + l, ABmz + k * l);
	//g = glm::vec3(ABmx + i * l + l, ABmy + j * l + l, ABmz + k * l + l);
	//h = glm::vec3(ABmx + i * l, ABmy + j * l + l, ABmz + k * l + l);

	ga = gradientFild[i * gj * gk + j * gk + k];
	gb = gradientFild[(i + 1) * gj * gk + j * gk + k];
	gc = gradientFild[(i + 1) * gj * gk + j * gk + k + 1];
	gd = gradientFild[i * gj * gk + j * gk + k + 1];

	ge = gradientFild[i * gj * gk + (j+1) * gk + k];
	gf = gradientFild[(i + 1) * gj * gk + (j + 1) * gk + k];
	gg = gradientFild[(i + 1) * gj * gk + (j + 1) * gk + k + 1];
	gh = gradientFild[i * gj * gk + (j + 1) * gk + k + 1];

	//interpolation in a b c d(i,k)
	float factori = (X[ind].x - a.x) / l;
	float factork = (X[ind].z - a.z) / l;
	float factorj = (X[ind].y - a.y) / l;

	glm::vec3 gabcd = (ga * (1 - factori) + gb * factori) * (1 - factork) + (gd * (1 - factori) + gc * factori) * factork;
	glm::vec3 gefgh = (ge * (1 - factori) + gf * factori) * (1 - factork) + (gh * (1 - factori) + gg * factori) * factork;

	*g = gabcd * (1 - factorj) + gefgh * factorj;

}

void DeformableObject::GenerateBlocks(size_t xdim, size_t ydim, size_t zdim, float width, float height, float depth) {
	total_points = (xdim + 1)*(ydim + 1)*(zdim + 1);
	X.resize(total_points);
	Xi.resize(total_points);
	IsFixed.resize(total_points);
	int ind = 0;
	float hzdim = zdim / 2.0f;
	for (size_t x = 0; x <= xdim; ++x) {
		for (unsigned int y = 0; y <= ydim; ++y) {
			for (unsigned int z = 0; z <= zdim; ++z) {
				X[ind] = glm::vec3(width*x, height*z + 0.5, depth*y);
				Xi[ind] = X[ind];

				//Make the first few points fixed
				//if (Xi[ind].x < 0.01)
					IsFixed[ind] = false;

				ind++;
			}
		}
	}
	//offset the tetrahedral mesh by 0.5 units on y axis
	//and 0.5 of the depth in z axis
	/*
	for (size_t i = 0; i<total_points; i++) {
		X[i].y += 0.5;
		X[i].z -= hzdim*depth;
	}
*/
	for (size_t i = 0; i < xdim; ++i) {
		for (size_t j = 0; j < ydim; ++j) {
			for (size_t k = 0; k < zdim; ++k) {
				int p0 = (i * (ydim + 1) + j) * (zdim + 1) + k;
				int p1 = p0 + 1;
				int p3 = ((i + 1) * (ydim + 1) + j) * (zdim + 1) + k;
				int p2 = p3 + 1;
				int p7 = ((i + 1) * (ydim + 1) + (j + 1)) * (zdim + 1) + k;
				int p6 = p7 + 1;
				int p4 = (i * (ydim + 1) + (j + 1)) * (zdim + 1) + k;
				int p5 = p4 + 1;
				// Ensure that neighboring tetras are sharing faces
				if ((i + j + k) % 2 == 1) {
					AddTetrahedron(p1, p2, p6, p3);
					AddTetrahedron(p3, p6, p4, p7);
					AddTetrahedron(p1, p4, p6, p5);
					AddTetrahedron(p1, p3, p4, p0);
					AddTetrahedron(p1, p6, p4, p3);
				}
				else {
					AddTetrahedron(p2, p0, p5, p1);
					AddTetrahedron(p2, p7, p0, p3);
					AddTetrahedron(p2, p5, p7, p6);
					AddTetrahedron(p0, p7, p5, p4);
					AddTetrahedron(p2, p0, p7, p5);
				}
				total_tetrahedra += 5;
			}
		}
	}

}
void DeformableObject::Reset()
{
	for (int i = 0; i < total_points; ++i)
	{
		X[i] = Xi[i];
	}
	memset(&(V[0].x), 0, total_points*sizeof(glm::vec3));
	debug << "Reset all." << endl;
	//debug2 << "Reset all." << endl;
}
void DeformableObject::firstPass(HashMap *H, int objId)
{
	for (int j = 0; j < total_points; ++j)
	{
		glm::vec3 p = X[j];
		int x = (int)(p.x / g_l);
		int y = (int)(p.y / g_l);
		int z = (int)(p.z / g_l);
		int h = ((x * 73856093) ^ (y * 19349663) ^ (z * 83492791)) % 199;
		if (h < 0) h = -h;

		if (H->cell[h].T != H->T)
		{
			H->cell[h].nodes.clear();
			H->cell[h].T = H->T;
		}
		Vertex v;
		v.objId = objId;
		v.localIndex = j;
		v.p = p;
		H->cell[h].nodes.push_back(v);
	}
}
void DeformableObject::secondPass(HashMap *H, int objId)
{//vector parameter  a vector is expensive;
	//debug << "---------------------" << 2 << endl;

	for (int i = 0; i < total_tetrahedra; ++i)
	{
		float Mx, My, Mz;//max
		float mx, my, mz;//min
		int n1 = tetrahedra[i].indices[0];
		int n2 = tetrahedra[i].indices[1];
		int n3 = tetrahedra[i].indices[2];
		int n4 = tetrahedra[i].indices[3];
		glm::vec3 v0 = X[n1];
		glm::vec3 v1 = X[n2];
		glm::vec3 v2 = X[n3];
		glm::vec3 v3 = X[n4];

		mx = glm::min(v0.x, v1.x);
		mx = glm::min(mx, v2.x);
		mx = glm::min(mx, v3.x);
		my = glm::min(v0.y, v1.y);
		my = glm::min(my, v2.y);
		my = glm::min(my, v3.y);
		mz = glm::min(v0.z, v1.z);
		mz = glm::min(mz, v2.z);
		mz = glm::min(mz, v3.z);
		Mx = glm::max(v0.x, v1.x);
		Mx = glm::max(Mx, v2.x);
		Mx = glm::max(Mx, v3.x);
		My = glm::max(v0.y, v1.y);
		My = glm::max(My, v2.y);
		My = glm::max(My, v3.y);
		Mz = glm::max(v0.z, v1.z);
		Mz = glm::max(Mz, v2.z);
		Mz = glm::max(Mz, v3.z);

        int M_x = Mx / g_l + 1;
		int M_y = My / g_l + 1;
		int M_z = Mz / g_l + 1;

		int m_x = mx / g_l;
		int m_y = my / g_l;
		int m_z = mz / g_l;

		
		//debug << mx << " " << my << " " << mz << " " << Mx << " " << My << " " << Mz << "              ";
		//debug << mx / g_l << " " << my / g_l << " " << mz / g_l << " " << Mx / g_l << " " << My / g_l << " " << Mz / g_l << "              ";
		//debug << m_x << " " << m_y << " " << m_z << " " << M_x << " " << M_y << " " << M_z << endl;

		for (int j = m_x; j <= M_x; ++j)
			for (int k = m_y; k <= M_y; ++k)
				for (int l = m_z; l <= M_z; ++l)
				{
			        int h = ((j * 73856093) ^ (k * 19349663) ^ (l * 83492791)) % 199;
			        if (h < 0) h = -h;
					
					//debug << h << " " << H->cell[h].T << "        " << H->T << endl;

			        if (H->cell[h].T == H->T)
					{
			        	list<Vertex> ::iterator node_iter;
			        	for (node_iter = H->cell[h].nodes.begin(); node_iter != H->cell[h].nodes.end(); node_iter++)
			        	{
			        		if (node_iter->objId == objId && ((node_iter->localIndex == n1) || (node_iter->localIndex == n2) || (node_iter->localIndex == n3) || (node_iter->localIndex == n4)))
			        			continue;
			        		glm::mat3 A;
			        		glm::vec3 e1, e2, e3;
			        		e1 = v1 - v0;
			        		e2 = v2 - v0;
			        		e3 = v3 - v0;
			        
			        		A[0][0] = e1.x;    A[0][1] = e2.x;    A[0][2] = e3.x;
			        		A[1][0] = e1.y;    A[1][1] = e2.y;    A[1][2] = e3.y;
			        		A[2][0] = e1.z;    A[2][1] = e2.z;    A[2][2] = e3.z;

			                A = glm::inverse(A);
			        		glm::vec3 beta = node_iter->p - v0;
							glm::vec3 b;
							//that is something wrong with beta = A * beta;
							b.x = A[0][0] * beta.x + A[0][1] * beta.y + A[0][2] * beta.z;
							b.y = A[1][0] * beta.x + A[1][1] * beta.y + A[1][2] * beta.z;
							b.z = A[2][0] * beta.x + A[2][1] * beta.y + A[2][2] * beta.z;

							//debug << b.x * 100 << "   " << b.y * 100 << "   " << b.z * 100 << endl;

			        		if (b.x >= 0 && b.y >= 0 && b.z >= 0 && (b.x + b.y + b.z) <= 1)
			        		{//contact(n,t)
			        			//glm::vec3 g0(0), g1(0), g2(0), g3(0);
			        
			        			float d = distanceV[n1] * (1 - b.x - b.y - b.z) + distanceV[n2] * b.x + distanceV[n3] * b.y + distanceV[n4] * b.z;
								
								glm::vec3 g = gradientV[n1] * (1 - beta.x - beta.y - beta.z) + gradientV[n2] * beta.x + gradientV[n3] * beta.y + gradientV[n4] * beta.z;
								g = tetrahedra[i].Re * g;
								//g = glm::vec3(0.0, 1.0f, 0.0f);
			        			glm::vec3 f = 250000.0f * d * g;

								F_add f_add;
								f_add.LocalId = node_iter->localIndex;
								f_add.f = f;
								g_F[node_iter->objId].push_back(f_add);


			        		}
			        
			        	}
			        }
				}
	}
}
void DeformableObject::F_ADD(int objId)
{
	for (int i = 0; i < g_F[objId].size(); ++i)
	{
		F[g_F[objId][i].LocalId] += g_F[objId][i].f;
	}
}

void DeformableObject::CalculateVN()
{
	for (int i = 0; i < bTriangle.size(); ++i)
	{
		int n0 = bTriangle[i].indices[0];
		int n1 = bTriangle[i].indices[1];
		int n2 = bTriangle[i].indices[2];
		glm::vec3 e01, e02;
		
		e01 = X[n1] - X[n0];
		e02 = X[n2] - X[n0];
		normaloftriangle[i] = glm::normalize(glm::cross(e02, e01));	
	}
	for (int i = 0; i < total_points; ++i)
	{
		glm::vec3 normal(0);
		for (int j = 0; j < trianglesOfVertices[i].size(); ++j)
		{
			normal += normaloftriangle[trianglesOfVertices[i][j]];
		}
		normalofvetex[i] = glm::normalize(normal);
	}
}
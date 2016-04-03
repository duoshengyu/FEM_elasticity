#include "deformableobject.h"
ofstream debug2("debug2.txt");
extern ofstream debug;
extern vector<DeformableObject> g_models;
extern float g_l;//global grid size
extern vector<vector<F_add>> g_F;
//-----------------------------------------------------------------------
//Constructe function. Include read mesh, initiate necessary buffer.
//Initiate FEM buffer, calculate distance fild and gradient fild
// G = (d(i+1,j,k)-d(i-1,j,k), d(i,j+1,k) - d(i,j-1,k), d(i,j,k+1) - d(i,j,k-1) )
//-----------------------------------------------------------------------
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
//-----------------------------------------------------------------------
// Read mesh from file
//-----------------------------------------------------------------------
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

//-----------------------------------------------------------------------
//Orthogonalization rotation matrix
//-----------------------------------------------------------------------

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

//-----------------------------------------------------------------------
//draw mesh
//-----------------------------------------------------------------------
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

//-----------------------------------------------------------------------
//Distance Fild: Shortest distance vertices in grid to mesh
//Gradient Fild:
// G = (d(i+1,j,k)-d(i-1,j,k), d(i,j+1,k) - d(i,j-1,k), d(i,j,k+1) - d(i,j,k-1) )
//-----------------------------------------------------------------------

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
//-----------------------------------------------------------------------
//To determine whether a vertex is in a tetrahedron
//If a vertex is in a tet then it can represent as the weight multiply
//sum of four vertices. s.t 0 < weight < 1.
//-----------------------------------------------------------------------
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
//-----------------------------------------------------------------------
//Shortest distant of a vertex to mesh
//-----------------------------------------------------------------------
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
//-----------------------------------------------------------------------
//Shortest distant of a vertex to a triangle
//-----------------------------------------------------------------------
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
//-----------------------------------------------------------------------
//calculte shortest distance and gradient of arbitrary point using interpolate
//-----------------------------------------------------------------------
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
//-----------------------------------------------------------------------
//Spatial hasing first pass 
//hasing all point to space, mark with mesh index
//-----------------------------------------------------------------------
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
//-----------------------------------------------------------------------
//Spatial hasing second pass 
//AABB a tetrahedron, search hash table and determine whether the vertex
//is in this tetrahedron
//-----------------------------------------------------------------------
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
		//AABB of a tethedron
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
		


		//search in hash table
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
//-----------------------------------------------------------------------
//Calculate vertex nomal
//-----------------------------------------------------------------------
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
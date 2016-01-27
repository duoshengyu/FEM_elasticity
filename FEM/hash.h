#pragma once
#include <App.h>
#include <list>
#include <vector>
struct Vertex
{
	int objId;
	int localIndex;
	glm::vec3 p;
};

class HashCell
{
public:
	double T = 0;
	list<Vertex> nodes;
	HashCell()
	{
	}
	~HashCell()
	{
	}
};

class HashMap
{
public:
	double T = 0;
	vector<HashCell> cell;

	HashMap(int n)
	{
		for (int i = 0; i < n; ++i)
		{
			HashCell hc;
			cell.push_back(hc);
		}
	}
	~HashMap()
	{
	}
};

#ifndef CAMERA_H
#define CAMERA_H
//--------------------------------------------------------------------------------------------
//simple track ball camera class for OpenGL
//--------------------------------------------------------------------------------------------
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
using namespace std;

class TrackballCamera
{
public:
	TrackballCamera(float w, float h);

	void                SetView(glm::vec3 pEye, glm::vec3 pLookat);
	void                SetProj(float fFov, float fAspect, float fNearPlane, float fFarPlane);
	void computeQuat();
	void computeTran();
	//get proj * view * world
	glm::mat4 getMVP()                                    { return mProj * mView * mmWorld; }
	glm::mat4 getM()                                      { return mmWorld; }
	glm::mat4 getV()                                      { return mView; }
	//calculate model matrix using parameters
	void setMmworldQuat()                                 { mmWorld = glm::mat4_cast(mQuat) * mmWorld; };
	void setMmworldTran()                                 { mmWorld = mTran * mmWorld; };
	void setMmworldScle()                                 { mmWorld = mScale * mmWorld; };
	//set scale factor
	void setScaleFactor(float x)                          { mSFactor = x; mScale = glm::scale(glm::mat4(),glm::vec3(x)); }
	void initMousePosition(float x, float y)              { SetCurMousePosition(x, y); SetPreMousePosition(x, y); }
	void SetMouseLButtonStat(bool stat)                   { mbMouseLButtonDown = stat; }
	void SetMouseLWheelStat(bool stat)                    { mbMouseWheelRoll = stat; }
	void SetMouseRButtonStat(bool stat)                   { mbMouseRButtonDown = stat; }

	void SetCurMousePosition(float x, float y)            { curMousePosition.x = x; curMousePosition.y = y; }
	void SetPreMousePosition(float x, float y)            { preMousePosition.x = x; preMousePosition.y = y; }

	void SetWindow(float w, float h)                      { windowWidth = w; windowHeight = h; }

	bool IsMouseLButtonDown() const                       { return mbMouseLButtonDown; }
	bool IsMouseWheelRoll() const                         { return mbMouseRButtonDown; }
	bool IsMouseRButtonDown() const                       { return mbMouseRButtonDown; }

protected:
	glm::vec3 mCenter;
	float mRadius;
	glm::quat mQuat;//rotate quaternion
	glm::mat4 mTran;//translate matrix
	glm::mat4 mScale;
	float     mSFactor;

	glm::mat4 mView;
	glm::mat4 mProj;
	glm::mat4 mmWorld;

	bool mbMouseLButtonDown;    // True if left button is down 
	bool mbMouseWheelRoll;          // True if middle wheel is roll 
	bool mbMouseRButtonDown;    // True if right button is down 

	float windowWidth;
	float windowHeight;
	glm::vec2 curMousePosition;
	glm::vec2 preMousePosition;
};

#endif 
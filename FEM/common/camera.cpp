#include "camera.h"

TrackballCamera::TrackballCamera(float w, float h)
{
	SetWindow(w, h);
	mCenter = glm::vec3(0.0f);
	float mRadius = 1.0f;
	mQuat = glm::quat();

	glm::vec3 pos = glm::vec3(0.0f, 0.0f, 1.0f);
	glm::vec3 target = glm::vec3(0.0f);

	//setMmworldQuat();
	mmWorld = glm::scale(glm::mat4(1.0f), glm::vec3(2.0f));
	glm::mat4 t = glm::translate(mTran, glm::vec3(0.0f, 0.0f, 0.0f));
	mmWorld *= t;
	SetView(pos, target);
	SetProj(45.0f, (float)windowWidth / windowHeight, 0.1f, 100.f);

	mbMouseLButtonDown = false;
	mbMouseWheelRoll = false;
	mbMouseRButtonDown = false;
}

void TrackballCamera::SetView(glm::vec3 eye, glm::vec3 target)
{
	glm::vec3 up(0.0f, 1.0f, 0.0f);
	mView = glm::lookAt(eye, target, up);
}

void TrackballCamera::SetProj(float fFov, float fAspect, float fNearPlane, float fFarPlane)
{
	mProj = glm::perspective(fFov, fAspect, fNearPlane, fFarPlane);
}
//calculate rotate quaternion
void TrackballCamera::computeQuat()
{
	float d, as;
	glm::vec3 preVec, curVec, vMove, axis;

	//transform from window coordinate to projection plane
	preVec.x = (2.0f * preMousePosition.x / windowWidth - 1.0f) / mProj[0][0];
	preVec.y = (-2.0f * preMousePosition.y / windowHeight + 1.0f) / mProj[1][1];
	preVec.z = 1.0f;
	d = preVec.x * preVec.x + preVec.y * preVec.y;

	curVec.x = (2.0f * curMousePosition.x / windowWidth - 1.0f) / mProj[0][0];
	curVec.y = (-2.0f * curMousePosition.y / windowHeight + 1.0f) / mProj[1][1];
	curVec.z = 1.0f;

	vMove = curVec - preVec;
	//calculate rotate axis
	axis = glm::cross(curVec, preVec);
	axis = glm::normalize(axis);

	as = (glm::length(vMove) * glm::length(vMove) - glm::length(curVec) * glm::length(curVec) - glm::length(preVec) * glm::length(preVec)) / (2.0f*glm::length(curVec)*glm::length(preVec));

	mQuat = glm::angleAxis(2.0f * as, axis);

	preMousePosition = curMousePosition;
}
void TrackballCamera::computeTran()
{
	glm::vec3 preVec, curVec, vMove;
	preVec.x = (2.0f * preMousePosition.x / windowWidth - 1.0f) / mProj[0][0];
	preVec.y = (-2.0f * preMousePosition.y / windowHeight + 1.0f) / mProj[1][1];
	preVec.z = 0.0f;

	curVec.x = (2.0f * curMousePosition.x / windowWidth - 1.0f) / mProj[0][0];
	curVec.y = (-2.0f * curMousePosition.y / windowHeight + 1.0f) / mProj[1][1];
	curVec.z = 0.0f;

	vMove = curVec - preVec;

	mTran[3][0] = 0.0f;
	mTran[3][1] = 0.0f;
	mTran[3][2] = 0.0f;
	mTran = glm::translate(mTran, vMove);

	preMousePosition = curMousePosition;
}
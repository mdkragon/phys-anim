#include "Tests.h"
#include <iostream>
#include <limits>
#include <string>
#include <gl/glut.h>
#include <OgreMath.h>
#include "World.h"
#include "RigidBody.h"
#include "Ground.h"
#include "Box.h"
#include "Sphere.h"
#define NOMINMAX
#include <windows.h>
#include <utility>
#include "Material.h"
#include "WorldLoader.h"
#include <windef.h>
#include <commdlg.h>
#include <MMSystem.h>
#include "SoundManager.h"

// fixes for weird compile error on illegal token on right side of :: for max and min
#undef min
#undef max

void SetupGraphics();
void RunGame();
void LaunchMissile();
void LoadWorldIntoGame( std::string const& mbfilename );

World g_world;


int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	std::cout.sync_with_stdio();

	//RunTests();
	SetupGraphics();

	// LOOK if you wanted, you could load a level right here for testing
	LoadWorldIntoGame("../Worlds/justground.xml");
	//LoadWorldIntoGame("../Worlds/seesaw.xml");
	//LoadWorldIntoGame("../Worlds/test_many.xml");

	g_world.SetUseSweepAndPrune(false);
	RunGame();
	
	/*
	g_world.SetUseSweepAndPrune(true);
	RunBenchmark("../Worlds/test_40x40.xml", 100);
	g_world.SetUseSweepAndPrune(false);
	RunBenchmark("../Worlds/test_40x40.xml", 100);
	*/
}

void display();
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void keyboard(unsigned char key, int x, int y);
void reshape(int width, int height);
void idle();

int g_width = 640;
int g_height = 480;
Vector3 g_cameraPos = Vector3(0,3,30);
// angles in radians
float g_cameraHeading = 0;
float g_cameraPitch = 0;
int g_activeButton = -1;
int g_lastX, g_lastY;
int g_prevTime;
Material* g_missileMaterial = NULL;
int g_stepNum = 0;

const float TRANSLATE_SPEED = 0.1f;
const float ROTATE_SPEED = 0.005f;

void SetupGraphics()
{
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(640, 480);
	glutCreateWindow("Boom Blox!");

	glutDisplayFunc(display);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutIdleFunc(idle);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1,1);
	glClearColor(0.4, 0.5, 1.0, 0.0);
}

void display()
{
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, float(g_width)/g_height, 0.1, 500.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(g_cameraPitch*180.0f/M_PI, -1, 0, 0);
	glRotatef(g_cameraHeading*180.0f/M_PI, 0, -1, 0);
	glTranslatef(-g_cameraPos.x, -g_cameraPos.y, -g_cameraPos.z);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	float lightPos[] = {-50, 40, 30, 0};
	float ambient[] = {0.4, 0.4, 0.4, 1};
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glEnable(GL_COLOR_MATERIAL);

	g_world.Render();

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, g_width, g_height, 0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glColor3f(1,1,1);
	glBegin(GL_POLYGON);
	glVertex2i(1,1);
	glVertex2i(16,1);
	glVertex2i(16,16);
	glVertex2i(1,16);
	glEnd();
	glColor3f(0,0,0);
	glBegin(GL_LINE_STRIP);
	glVertex2i(5, 3);
	glVertex2i(5, 14);
	glVertex2i(11, 14);
	glEnd();
	glutSwapBuffers();
}

void LoadWorldIntoGame( std::string const& mbfilename )
{
	g_world = LoadWorldFromFile(mbfilename);
	g_missileMaterial = new Material(Material::MISSILE_STUFF);
	g_world.AddMaterial(g_missileMaterial);
}

void PromptToLoad()
{
	OPENFILENAME ofn;
	wchar_t filename[1024];
	memset(filename, 0, sizeof(filename));

	memset(&ofn, 0, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.lpstrFilter = L"World Files\0*.xml\0All Files\0*.*\0\0";
	ofn.lpstrFile = filename;
	ofn.nMaxFile = 1024;
	ofn.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY;

	BOOL result = GetOpenFileName(&ofn);
	if(result)
	{
		char mbfilename[1024];
		wcstombs(mbfilename, filename, 1024);
		LoadWorldIntoGame(mbfilename);

	}
	else
	{
		DWORD err = CommDlgExtendedError();
		err = err;
	}

	g_stepNum = 0;
}
void mouse(int button, int state, int x, int y)
{
	if(state == GLUT_DOWN)
	{
		if(x < 16 && y < 16)
		{
			PromptToLoad();
		}
		else
		{
			g_activeButton = button;
			g_lastX = x;
			g_lastY = y;
		}
	}
	else if(state == GLUT_UP)
	{
		g_activeButton = -1;
	}
}

void motion(int x, int y)
{
	switch(g_activeButton)
	{
	case GLUT_LEFT_BUTTON:
		{
			Vector3 trans((x-g_lastX)*TRANSLATE_SPEED, 0, (y-g_lastY)*TRANSLATE_SPEED);
			trans = Quaternion(g_cameraHeading, Vector3::UNIT_Y) * trans;
			g_cameraPos += trans;
			break;
		}

	case GLUT_RIGHT_BUTTON:
		{
			g_cameraHeading -= (x-g_lastX)*ROTATE_SPEED;
			g_cameraPitch -= (y-g_lastY)*ROTATE_SPEED;
			break;
		}
	default:
		return;
	}
	g_lastX = x;
	g_lastY = y;
}

Vector3 ComputeTrajectory()
{
	// LOOK for better missile control
	const float MISSILE_SPEED = 40;
	float heading = g_cameraHeading;
	float pitch = g_cameraPitch;

	return Vector3(
		-sinf(heading)*cosf(pitch),
		sinf(pitch),
		-cosf(heading)*cosf(pitch)) * MISSILE_SPEED;
}

void LaunchMissile()
{
	Sphere* s = new Sphere(1);
	s->SetMaterial(g_missileMaterial);
	s->SetPosition(g_cameraPos);
	s->SetVelocity(ComputeTrajectory());

	g_world.AddBody(s);
}

void keyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
	case ' ':
		LaunchMissile();
		break;
	}
}

void reshape(int width, int height)
{
	g_width = width;
	g_height = height;
	
	glViewport(0, 0, width, height);
}

void idle()
{
	int newTime = timeGetTime();
	float dT = (newTime-g_prevTime) / 1000.0f;
	g_prevTime = newTime;

	dT = std::min(dT, 0.03f);

	
	// LOOK fixed timesteps are much easier to debug
	g_world.Simulate(dT);
	//g_world.Simulate(0.03f);

	g_stepNum++;
	glutPostRedisplay();
}

void RunGame()
{
	timeBeginPeriod(1);
	g_prevTime = timeGetTime();

	glutMainLoop();
}



















///*===============================================================================================
// 3d Example
// Copyright (c), Firelight Technologies Pty, Ltd 2004-2011.
//
// This example shows how to basic 3d positioning
//===============================================================================================*/
//
//#include <windows.h>
//#include <stdio.h>
//#include <conio.h>
//#include <math.h>
//
//#include "fmod.hpp"
//#include "fmod_errors.h"
//
//const int   INTERFACE_UPDATETIME = 50;      // 50ms update for interface
//const float DISTANCEFACTOR = 1.0f;          // Units per meter.  I.e feet would = 3.28.  centimeters would = 100.
//
//
//void ERRCHECK(FMOD_RESULT result)
//{
//    if (result != FMOD_OK)
//    {
//        printf("FMOD error! (%d) %s\n", result, FMOD_ErrorString(result));
//        exit(-1);
//    }
//}
//
//
//int main(int argc, char *argv[])
//{
//    FMOD::System    *system;
//    FMOD::Sound     *sound1, *sound2, *sound3;
//    FMOD::Channel   *channel1 = 0, *channel2 = 0, *channel3 = 0;
//    FMOD_RESULT      result;
//    int              key, numdrivers;
//    bool             listenerflag = true;
//    FMOD_VECTOR      listenerpos  = { 0.0f, 0.0f, -1.0f * DISTANCEFACTOR };
//    unsigned int     version;
//    FMOD_SPEAKERMODE speakermode;
//    FMOD_CAPS        caps;
//    char             name[256];
//
//    printf("===============================================================\n");
//    printf("3d Example.  Copyright (c) Firelight Technologies 2004-2011.\n");
//    printf("===============================================================\n");
//    printf("This example plays 2 3D sounds in hardware.  Optionally you can\n");
//    printf("play a 2D hardware sound as well.\n");
//    printf("===============================================================\n\n");
//
//    /*
//        Create a System object and initialize.
//    */
//    result = FMOD::System_Create(&system);
//    ERRCHECK(result);
//    
//    result = system->getVersion(&version);
//    ERRCHECK(result);
//
//    if (version < FMOD_VERSION)
//    {
//        printf("Error!  You are using an old version of FMOD %08x.  This program requires %08x\n", version, FMOD_VERSION);
//        return 0;
//    }
//    
//    result = system->getNumDrivers(&numdrivers);
//    ERRCHECK(result);
//
//    if (numdrivers == 0)
//    {
//        result = system->setOutput(FMOD_OUTPUTTYPE_NOSOUND);
//        ERRCHECK(result);
//    }
//    else
//    {
//        result = system->getDriverCaps(0, &caps, 0, &speakermode);
//        ERRCHECK(result);
//
//        result = system->setSpeakerMode(speakermode);       /* Set the user selected speaker mode. */
//        ERRCHECK(result);
//
//        if (caps & FMOD_CAPS_HARDWARE_EMULATED)             /* The user has the 'Acceleration' slider set to off!  This is really bad for latency!. */
//        {                                                   /* You might want to warn the user about this. */
//            result = system->setDSPBufferSize(1024, 10);
//            ERRCHECK(result);
//        }
//
//        result = system->getDriverInfo(0, name, 256, 0);
//        ERRCHECK(result);
//
//        if (strstr(name, "SigmaTel"))   /* Sigmatel sound devices crackle for some reason if the format is PCM 16bit.  PCM floating point output seems to solve it. */
//        {
//            result = system->setSoftwareFormat(48000, FMOD_SOUND_FORMAT_PCMFLOAT, 0,0, FMOD_DSP_RESAMPLER_LINEAR);
//            ERRCHECK(result);
//        }
//    }
//
//    result = system->init(100, FMOD_INIT_NORMAL, 0);
//    if (result == FMOD_ERR_OUTPUT_CREATEBUFFER)         /* Ok, the speaker mode selected isn't supported by this soundcard.  Switch it back to stereo... */
//    {
//        result = system->setSpeakerMode(FMOD_SPEAKERMODE_STEREO);
//        ERRCHECK(result);
//            
//        result = system->init(100, FMOD_INIT_NORMAL, 0);/* ... and re-init. */
//        ERRCHECK(result);
//    }
//
//    
//    /*
//        Set the distance units. (meters/feet etc).
//    */
//    result = system->set3DSettings(1.0, DISTANCEFACTOR, 1.0f);
//    ERRCHECK(result);
//
//    /*
//        Load some sounds
//    */
//    result = system->createSound("drumloop.wav", FMOD_3D, 0, &sound1);
//    ERRCHECK(result);
//    result = sound1->set3DMinMaxDistance(0.5f * DISTANCEFACTOR, 5000.0f * DISTANCEFACTOR);
//    ERRCHECK(result);
//    result = sound1->setMode(FMOD_LOOP_NORMAL);
//    ERRCHECK(result);
//
//    result = system->createSound("jaguar.wav", FMOD_3D, 0, &sound2);
//    ERRCHECK(result);
//    result = sound2->set3DMinMaxDistance(0.5f * DISTANCEFACTOR, 5000.0f * DISTANCEFACTOR);
//    ERRCHECK(result);
//    result = sound2->setMode(FMOD_LOOP_NORMAL);
//    ERRCHECK(result);
//
//    result = system->createSound("swish.wav", FMOD_SOFTWARE | FMOD_2D, 0, &sound3);
//    ERRCHECK(result);
//
//    /*
//        Play sounds at certain positions
//    */
//    {
//        FMOD_VECTOR pos = { -10.0f * DISTANCEFACTOR, 0.0f, 0.0f };
//        FMOD_VECTOR vel = {  0.0f, 0.0f, 0.0f };
//
//        result = system->playSound(FMOD_CHANNEL_FREE, sound1, true, &channel1);
//        ERRCHECK(result);
//        result = channel1->set3DAttributes(&pos, &vel);
//        ERRCHECK(result);
//        result = channel1->setPaused(false);
//        ERRCHECK(result);
//    }
//
//    {
//        FMOD_VECTOR pos = { 15.0f * DISTANCEFACTOR, 0.0f, 0.0f };
//        FMOD_VECTOR vel = { 0.0f, 0.0f, 0.0f };
//
//        result = system->playSound(FMOD_CHANNEL_FREE, sound2, true, &channel2);
//        ERRCHECK(result);
//        result = channel2->set3DAttributes(&pos, &vel);
//        ERRCHECK(result);
//        result = channel2->setPaused(false);
//        ERRCHECK(result);
//    }
//
//    /*
//        Display help
//    */
//    {
//        int numchannels;
//
//        result = system->getHardwareChannels(&numchannels);
//        ERRCHECK(result);
//    
//        printf("Hardware channels : %d\n", numchannels);
//    }
//
//    printf("=========================================================================\n");
//    printf("Press 1        Pause/Unpause 16bit 3D sound at any time\n");
//    printf("      2        Pause/Unpause 8bit 3D sound at any time\n");
//    printf("      3        Play 16bit STEREO 2D sound at any time\n");
//    printf("      <        Move listener left (in still mode)\n");
//    printf("      >        Move listener right (in still mode)\n");
//    printf("      SPACE    Stop/Start listener automatic movement\n");
//    printf("      ESC      Quit\n");
//    printf("=========================================================================\n");
//
//    /*
//        Main loop
//    */
//    do
//    {
//        if (_kbhit())
//        {
//            key = _getch();
//
//            if (key == '1') 
//            {
//                bool paused;
//                channel1->getPaused(&paused);
//                channel1->setPaused(!paused);
//            }
//
//            if (key == '2') 
//            {
//                bool paused;
//                channel2->getPaused(&paused);
//                channel2->setPaused(!paused);
//            }
//
//            if (key == '3') 
//            {
//                result = system->playSound(FMOD_CHANNEL_FREE, sound3, false, &channel3);
//                ERRCHECK(result);
//            }
//
//            if (key == ' ')
//            {
//                listenerflag = !listenerflag;
//            }
//
//            if (!listenerflag)
//            {
//                if (key == '<') 
//                {
//                    listenerpos.x -= 1.0f * DISTANCEFACTOR;
//                    if (listenerpos.x < -35 * DISTANCEFACTOR)
//                    {
//                        listenerpos.x = -35 * DISTANCEFACTOR;
//                    }
//                }
//                if (key == '>') 
//                {
//                    listenerpos.x += 1.0f * DISTANCEFACTOR;
//                    if (listenerpos.x > 36 * DISTANCEFACTOR)
//                    {
//                        listenerpos.x = 36 * DISTANCEFACTOR;
//                    }
//                }
//            }
//        }
//
//        // ==========================================================================================
//        // UPDATE THE LISTENER
//        // ==========================================================================================
//        {
//            static float t = 0;
//            static FMOD_VECTOR lastpos = { 0.0f, 0.0f, 0.0f };
//            FMOD_VECTOR forward        = { 0.0f, 0.0f, 1.0f };
//            FMOD_VECTOR up             = { 0.0f, 1.0f, 0.0f };
//            FMOD_VECTOR vel;
//
//            if (listenerflag)
//            {
//                listenerpos.x = (float)sin(t * 0.05f) * 33.0f * DISTANCEFACTOR; // left right pingpong
//            }
//
//            // ********* NOTE ******* READ NEXT COMMENT!!!!!
//            // vel = how far we moved last FRAME (m/f), then time compensate it to SECONDS (m/s).
//            vel.x = (listenerpos.x - lastpos.x) * (1000 / INTERFACE_UPDATETIME);
//            vel.y = (listenerpos.y - lastpos.y) * (1000 / INTERFACE_UPDATETIME);
//            vel.z = (listenerpos.z - lastpos.z) * (1000 / INTERFACE_UPDATETIME);
//
//            // store pos for next time
//            lastpos = listenerpos;
//
//            result = system->set3DListenerAttributes(0, &listenerpos, &vel, &forward, &up);
//            ERRCHECK(result);
//
//            t += (30 * (1.0f / (float)INTERFACE_UPDATETIME));    // t is just a time value .. it increments in 30m/s steps in this example
//
//            // print out a small visual display
//            {
//                char s[80];
//
//                sprintf(s, "|.......................<1>......................<2>....................|");
//
//                s[(int)(listenerpos.x / DISTANCEFACTOR) + 35] = 'L';
//                printf("%s\r", s);
//            }
//        }
//
//        system->update();
//
//        Sleep(INTERFACE_UPDATETIME - 1);
//    } while (key != 27);
//
//    printf("\n");
//
//    /*
//        Shut down
//    */
//    result = sound1->release();
//    ERRCHECK(result);
//    result = sound2->release();
//    ERRCHECK(result);
//    result = sound3->release();
//    ERRCHECK(result);
//
//    result = system->close();
//    ERRCHECK(result);
//    result = system->release();
//    ERRCHECK(result);
//
//    return 0;
//}
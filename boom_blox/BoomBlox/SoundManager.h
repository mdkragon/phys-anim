#ifndef SOUND_MANAGER_H
#define SOUND_MANAGER_H

#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <vector>
#include <string>
#include <Eigen\Dense>

#include "fmod.hpp"
#include "fmod_errors.h"

class SoundManager
{
public:
	SoundManager();
	~SoundManager();
	
	void Update();

	void PlayTestSound();
private:
	FMOD::System    *system;
    FMOD::Sound     *sound1, *sound2, *sound3;
    FMOD::Channel   *channel1, *channel2, *channel3;
    FMOD_SPEAKERMODE speakermode;
    FMOD_CAPS        caps;
    unsigned int     version;
    char             name[256];
    int              key, numdrivers;
    bool             listenerflag;// = true;
    FMOD_VECTOR      listenerpos;//  = { 0.0f, 0.0f, -1.0f * DISTANCEFACTOR };
    

};

#endif
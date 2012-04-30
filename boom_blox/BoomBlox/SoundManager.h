#ifndef SOUND_MANAGER_H
#define SOUND_MANAGER_H

#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <vector>
#include <string>
#include <Eigen\Dense>

//#include "Global.h"
#include "fmod.hpp"
#include "fmod_errors.h"

class SoundManager
{
public:
	SoundManager();
	~SoundManager();

	void Update();

	void PlayTestSound();
	void InitUserCreatedSound();
	void PlayUserCreatedSound();

	void InitCollisionSound();

  //FMOD_RESULT F_CALLBACK pcmsetposcallback(FMOD_SOUND *sound, int subsound, unsigned int position, FMOD_TIMEUNIT postype);

private:
	FMOD::System    *system;
    FMOD::Sound     *sound1, *sound2, *sound3, *usersound;
    FMOD::Channel   *channel1, *channel2, *channel3, *channeluser;
    FMOD_SPEAKERMODE speakermode;
    FMOD_CAPS        caps;


	// vector of sounds;
	std::vector<FMOD::Sound *> *sounds;
	std::vector<FMOD::Channel *> *channels;

    unsigned int     version;
    char             name[256];
    int              key, numdrivers;
    bool             listenerflag;// = true;
    FMOD_VECTOR      listenerpos;//  = { 0.0f, 0.0f, -1.0f * DISTANCEFACTOR };
	
};


FMOD_RESULT F_CALLBACK pcmsetposcallback(FMOD_SOUND *sound, int subsound, unsigned int position, FMOD_TIMEUNIT postype);
FMOD_RESULT F_CALLBACK pcmreadcallback(FMOD_SOUND *sound, void *data, unsigned int datalen);

//extern SoundManager* Sound_Manager;

#endif

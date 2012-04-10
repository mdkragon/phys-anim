#include "SoundManager.h"

const int   INTERFACE_UPDATETIME = 50;      // 50ms update for interface
const float DISTANCEFACTOR = 1.0f;          // Units per meter.  I.e feet would = 3.28.  centimeters would = 100.


void ERRCHECK(FMOD_RESULT result)
{
    if (result != FMOD_OK)
    {
        printf("FMOD error! (%d) %s\n", result, FMOD_ErrorString(result));
        exit(-1);
    }
}

SoundManager::SoundManager()
{
	/*
        Create a System object and initialize.
    */
	FMOD_RESULT result;
    result = FMOD::System_Create(&system);
    ERRCHECK(result);
    
    result = system->getVersion(&version);
    ERRCHECK(result);

    if (version < FMOD_VERSION)
    {
        printf("Error!  You are using an old version of FMOD %08x.  This program requires %08x\n", version, FMOD_VERSION);
        exit(1);
    }
    
    result = system->getNumDrivers(&numdrivers);
    ERRCHECK(result);

    if (numdrivers == 0)
    {
        result = system->setOutput(FMOD_OUTPUTTYPE_NOSOUND);
        ERRCHECK(result);
    }
    else
    {
        result = system->getDriverCaps(0, &caps, 0, &speakermode);
        ERRCHECK(result);

        result = system->setSpeakerMode(speakermode);       /* Set the user selected speaker mode. */
        ERRCHECK(result);

        if (caps & FMOD_CAPS_HARDWARE_EMULATED)             /* The user has the 'Acceleration' slider set to off!  This is really bad for latency!. */
        {                                                   /* You might want to warn the user about this. */
            result = system->setDSPBufferSize(1024, 10);
            ERRCHECK(result);
        }

        result = system->getDriverInfo(0, name, 256, 0);
        ERRCHECK(result);

        if (strstr(name, "SigmaTel"))   /* Sigmatel sound devices crackle for some reason if the format is PCM 16bit.  PCM floating point output seems to solve it. */
        {
            result = system->setSoftwareFormat(48000, FMOD_SOUND_FORMAT_PCMFLOAT, 0,0, FMOD_DSP_RESAMPLER_LINEAR);
            ERRCHECK(result);
        }
    }

    result = system->init(100, FMOD_INIT_NORMAL, 0);
    if (result == FMOD_ERR_OUTPUT_CREATEBUFFER)
    {
		/* Ok, the speaker mode selected isn't supported by this soundcard.  Switch it back to stereo... */
        result = system->setSpeakerMode(FMOD_SPEAKERMODE_STEREO);
        ERRCHECK(result);
            
		/* ... and re-init. */
        result = system->init(100, FMOD_INIT_NORMAL, 0);
        ERRCHECK(result);
    }

	/*
        Set the distance units. (meters/feet etc).
    */
    result = system->set3DSettings(1.0, DISTANCEFACTOR, 1.0f);
    ERRCHECK(result);


	/*
        Load some sounds
    */
    result = system->createSound("drumloop.wav", FMOD_3D, 0, &sound1);
    ERRCHECK(result);
    result = sound1->set3DMinMaxDistance(0.5f * DISTANCEFACTOR, 5000.0f * DISTANCEFACTOR);
    ERRCHECK(result);
    result = sound1->setMode(FMOD_LOOP_NORMAL);
    ERRCHECK(result);

    result = system->createSound("jaguar.wav", FMOD_3D, 0, &sound2);
    ERRCHECK(result);
    result = sound2->set3DMinMaxDistance(0.5f * DISTANCEFACTOR, 5000.0f * DISTANCEFACTOR);
    ERRCHECK(result);
    result = sound2->setMode(FMOD_LOOP_NORMAL);
    ERRCHECK(result);

    result = system->createSound("Bounce.wav", FMOD_SOFTWARE | FMOD_2D, 0, &sound3);
    ERRCHECK(result);
	
	{
        FMOD_VECTOR pos = { -10.0f * DISTANCEFACTOR, 0.0f, 0.0f };
        FMOD_VECTOR vel = {  0.0f, 0.0f, 0.0f };

        result = system->playSound(FMOD_CHANNEL_FREE, sound1, true, &channel1);
        ERRCHECK(result);
        result = channel1->set3DAttributes(&pos, &vel);
        ERRCHECK(result);
        result = channel1->setPaused(false);
        ERRCHECK(result);
    }
	
}

SoundManager::~SoundManager()
{
	/*
        Shut down
    */
	FMOD_RESULT result;
    result = sound1->release();
    ERRCHECK(result);
    result = sound2->release();
    ERRCHECK(result);
    result = sound3->release();
    ERRCHECK(result);

    result = system->close();
    ERRCHECK(result);
    result = system->release();
    ERRCHECK(result);
}

void SoundManager::Update()
{
	system->update();
}

void SoundManager::PlayTestSound()
{
	FMOD_RESULT result;

	result = system->playSound(FMOD_CHANNEL_FREE, sound3, false, &channel3);
	ERRCHECK(result);
	
	
 //   /*
 //       Play sounds at certain positions
 //   */
 //   {
 //       FMOD_VECTOR pos = { -10.0f * DISTANCEFACTOR, 0.0f, 0.0f };
 //       FMOD_VECTOR vel = {  0.0f, 0.0f, 0.0f };

 //       result = system->playSound(FMOD_CHANNEL_FREE, sound1, true, &channel1);
 //       ERRCHECK(result);
 //       result = channel1->set3DAttributes(&pos, &vel);
 //       ERRCHECK(result);
 //       result = channel1->setPaused(false);
 //       ERRCHECK(result);
 //   }

 //   {
 //       FMOD_VECTOR pos = { 15.0f * DISTANCEFACTOR, 0.0f, 0.0f };
 //       FMOD_VECTOR vel = { 0.0f, 0.0f, 0.0f };

 //       result = system->playSound(FMOD_CHANNEL_FREE, sound2, true, &channel2);
 //       ERRCHECK(result);
 //       result = channel2->set3DAttributes(&pos, &vel);
 //       ERRCHECK(result);
 //       result = channel2->setPaused(false);
 //       ERRCHECK(result);
 //   }
}
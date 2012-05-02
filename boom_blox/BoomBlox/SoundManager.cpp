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

SoundManager::SoundManager(int userCreatedFrequency)
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

    result = system->createSound("Bounce.wav", FMOD_SOFTWARE | FMOD_3D, 0, &sound3);
    ERRCHECK(result);

	m_userCreatedFrequency = userCreatedFrequency;
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

void SoundManager::SetListenerPose(Vector3 pos, float pitch, float heading)
{
	// only set attributes if they have changed
	if (pos == m_listenerPosition && pitch == m_listenerPitch && heading == m_listenerHeading) {
		return;
	}

	// TODO: is this correct...
	//  graphics axis is: x-right, y-up, z-out of screen
	//  sound axis uses left handed: x-right, y-up, z-in to screen
	//    forward vector is -z axis of camera
	//    up vector is +y of camera
	
	// create rotation matrix
	Matrix3 Rpitch;
	Matrix3 Rheading;
	Rpitch.FromAxisAngle(Vector3(1.0f, 0.0f, 0.0f), pitch);
	Rheading.FromAxisAngle(Vector3(0.0f, 1.0f, 0.0f), heading);
	Vector3 forward = Rpitch * Rheading * Vector3(0.0f, 0.0f, -1.0f);
	Vector3 up = Rpitch * Rheading * Vector3(0.0f, 1.0f, 0.0f);

	
	// negate z component for left handed coordinate system
	FMOD_VECTOR listenerPosition = {pos.x, pos.y, -pos.z};
	FMOD_VECTOR listenerForward = {forward.x, forward.y, -forward.z};
	FMOD_VECTOR listenerUp = {up.x, up.y, -up.z};

	/*
	printf("setting listener pose to: (%04.2f, %04.2f, %04.2f), (%04.2f, %04.2f, %04.2f), (%04.2f, %04.2f, %04.2f)\n",
				listenerPosition.x, listenerPosition.y, listenerPosition.z, 
				listenerForward.x, listenerForward.y, listenerForward.z,
				listenerUp.x, listenerUp.y, listenerUp.z); 
    */

	FMOD_RESULT result;
	result = system->set3DListenerAttributes(0, &listenerPosition, 0, &listenerForward, &listenerUp);
	ERRCHECK(result);

	// store current pose
	m_listenerPosition = pos;
	m_listenerPitch = pitch;
	m_listenerHeading = heading;
}

void SoundManager::Update()
{
	system->update();
}

void SoundManager::PlayTestSound(Vector3 pos, Vector3 vel)
{
	FMOD_RESULT result;

	// convert position and velocity to FMOD vectors
	//	fmod uses left handed coordinate system so negate z components
	FMOD_VECTOR fpos = {pos.x, pos.y, -pos.z};
	FMOD_VECTOR fvel = {vel.x, vel.y, -vel.z};

	result = system->playSound(FMOD_CHANNEL_FREE, sound3, true, &channel3);
	ERRCHECK(result);
	// set the 3D position
    result = channel3->set3DAttributes(&fpos, &fvel);
    ERRCHECK(result);
	// unpause (play) the sound
    result = channel3->setPaused(false);
	ERRCHECK(result);
}

int SoundManager::GetUserCreatedFrequency()
{
	return m_userCreatedFrequency;
}

void SoundManager::PlayUserCreatedSample(float *data, int length, Vector3 pos, Vector3 vel)
{
	FMOD_RESULT result;
	FMOD::Sound *sound;
	FMOD::Channel *channel;
	int num_channels = 1;

	// convert position and velocity to FMOD vectors
	//	fmod uses left handed coordinate system so negate z components
	FMOD_VECTOR fpos = {pos.x, pos.y, -pos.z};
	FMOD_VECTOR fvel = {vel.x, vel.y, -vel.z};

	// sound info struct
	FMOD_CREATESOUNDEXINFO  createsoundexinfo;
	// sound mode
	FMOD_MODE mode = FMOD_3D | FMOD_OPENUSER | FMOD_LOOP_OFF | FMOD_SOFTWARE | FMOD_CREATESAMPLE | FMOD_OPENRAW | FMOD_OPENMEMORY;
	
	// fill the sound info struct
	memset(&createsoundexinfo, 0, sizeof(FMOD_CREATESOUNDEXINFO));
	// required
    createsoundexinfo.cbsize            = sizeof(FMOD_CREATESOUNDEXINFO);
	// Chunk size of stream update in samples.  This will be the amount of data passed to the user callback.
	createsoundexinfo.decodebuffersize  = m_userCreatedFrequency;
	// Length of PCM data in bytes of whole song (for Sound::getLength)
	createsoundexinfo.length            = length * sizeof(float);
	// Number of channels in the sound.
    createsoundexinfo.numchannels       = num_channels;
	// Default playback rate of sound.
    createsoundexinfo.defaultfrequency  = m_userCreatedFrequency;
	// Data format of sound.
    createsoundexinfo.format            = FMOD_SOUND_FORMAT_PCMFLOAT;
	// User callback for reading (we do not use this)
    createsoundexinfo.pcmreadcallback  	= NULL;
	// User callback for seeking (we do not use this)
    createsoundexinfo.pcmsetposcallback = NULL;

	// create the sound
	//   FMOD_RESULT F_API createSound (const char *name_or_data, FMOD_MODE mode, FMOD_CREATESOUNDEXINFO *exinfo, Sound **sound);
  	result = system->createSound((const char *)data, mode, &createsoundexinfo, &sound);
	ERRCHECK(result);


	// start the sound paused
	result = system->playSound(FMOD_CHANNEL_FREE, sound, true, &channel);
	ERRCHECK(result);
	// set the 3D position
    result = channel->set3DAttributes(&fpos, &fvel);
    ERRCHECK(result);
	// unpause (play) the sound
    result = channel->setPaused(false);
	ERRCHECK(result);
}


void SoundManager::fileLoader(){
       std::string filename = "matlab/rawr.txt";
       float curr_freq;
       std::filebuf fb;

       fb.open(filename, std::ios::in);
       std::istream is(&fb);

       for (int i = 0; i < 22050; i++) {
               is >> curr_freq;
               stream[i] = curr_freq;
       }
       fb.close();
}
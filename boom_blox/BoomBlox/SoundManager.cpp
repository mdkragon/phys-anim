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
	
	/*
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
	*/
	

	sounds = new std::vector<FMOD::Sound *>;
	channels = new std::vector<FMOD::Channel *>;
	
//	InitUserCreatedSound();

	sounds = new std::vector<FMOD::Sound *>;
    channels = new std::vector<FMOD::Channel *>;

	int length = 22050;
    stream = (float*)malloc(sizeof (float) * length);

	/*
    fileLoader();

	InitUserCreatedSample(stream, length);
	*/
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

void SoundManager::InitUserCreatedSound()
{
	// sound info struct
	FMOD_CREATESOUNDEXINFO  createsoundexinfo;
	// sound mode
	//FMOD_MODE mode = FMOD_3D | FMOD_OPENUSER | FMOD_LOOP_OFF | FMOD_HARDWARE;
	FMOD_MODE mode = FMOD_3D | FMOD_OPENUSER | FMOD_LOOP_OFF | FMOD_SOFTWARE;
	// actual sound object
	//FMOD::Sound *sound;
	int num_channels = 2;
	FMOD_RESULT result;

	FMOD::Sound *tmpsound;
	FMOD::Channel *tmpchannel;

	// fill the sound info struct
	memset(&createsoundexinfo, 0, sizeof(FMOD_CREATESOUNDEXINFO));
	// required
    createsoundexinfo.cbsize            = sizeof(FMOD_CREATESOUNDEXINFO);
	// Chunk size of stream update in samples.  This will be the amount of data passed to the user callback.
	createsoundexinfo.decodebuffersize  = 44100;
	// Length of PCM data in bytes of whole song (for Sound::getLength)
    createsoundexinfo.length            = 44100 * num_channels * sizeof(signed short) * 5;
	// Number of channels in the sound.
    createsoundexinfo.numchannels       = num_channels;
	// Default playback rate of sound.
    createsoundexinfo.defaultfrequency  = 44100;
	// Data format of sound.
    createsoundexinfo.format            = FMOD_SOUND_FORMAT_PCM16;
	// User callback for reading
    createsoundexinfo.pcmreadcallback  	= pcmreadcallback;
	// User callback for seeking.
    createsoundexinfo.pcmsetposcallback = pcmsetposcallback;

	// create the sound
	//   FMOD_RESULT F_API createSound (const char *name_or_data, FMOD_MODE mode, FMOD_CREATESOUNDEXINFO *exinfo, Sound **sound);
    result = system->createSound(0, mode, &createsoundexinfo, &tmpsound);
    ERRCHECK(result);
	sounds->push_back(tmpsound);
	channels->push_back(tmpchannel);
}

void SoundManager::InitUserCreatedSample(float *data, int length)
{
	// sound info struct
	FMOD_CREATESOUNDEXINFO  createsoundexinfo;
	// sound mode
	FMOD_MODE mode = FMOD_3D | FMOD_OPENUSER | FMOD_LOOP_OFF | FMOD_SOFTWARE | FMOD_CREATESAMPLE | FMOD_OPENRAW | FMOD_OPENMEMORY;
	

	int num_channels = 1;
	FMOD_RESULT result;

	FMOD::Sound *tmpsound;

	// fill the sound info struct
	memset(&createsoundexinfo, 0, sizeof(FMOD_CREATESOUNDEXINFO));
	// required
    createsoundexinfo.cbsize            = sizeof(FMOD_CREATESOUNDEXINFO);
	// Chunk size of stream update in samples.  This will be the amount of data passed to the user callback.
	createsoundexinfo.decodebuffersize  = 44100;
	// Length of PCM data in bytes of whole song (for Sound::getLength)
	createsoundexinfo.length            = length * sizeof(float);
	// Number of channels in the sound.
    createsoundexinfo.numchannels       = num_channels;
	// Default playback rate of sound.
    createsoundexinfo.defaultfrequency  = 44100;
	// Data format of sound.
    createsoundexinfo.format            = FMOD_SOUND_FORMAT_PCMFLOAT;
	// User callback for reading
    createsoundexinfo.pcmreadcallback  	= NULL;
	// User callback for seeking.
    createsoundexinfo.pcmsetposcallback = NULL;

	
	signed short * sdata = (signed short *)malloc(length*sizeof(signed short));
	for (int i = 0; i < length; i++) {
		sdata[i] = (signed short)(255 * data[i]);
	}
	

	printf("creating sound from data: %p\n", data);

	// create the sound
	//   FMOD_RESULT F_API createSound (const char *name_or_data, FMOD_MODE mode, FMOD_CREATESOUNDEXINFO *exinfo, Sound **sound);
    //result = system->createSound((const char *)sdata, mode, &createsoundexinfo, &tmpsound);
	result = system->createSound((const char *)data, mode, &createsoundexinfo, &tmpsound);
	ERRCHECK(result);


	// play the sound
	printf("playing sound\n");
	result = system->playSound(FMOD_CHANNEL_FREE, tmpsound, false, NULL);

	//Sleep(2);
}

void SoundManager::PlayUserCreatedSound() {
	FMOD_RESULT result;

	while(!(sounds->empty())) { // while sounds vector isn't empty
		FMOD::Sound * current = sounds->at(sounds->size()-1);
		FMOD::Channel * curChannel = channels->at(channels->size()-1);
		result = system->playSound(FMOD_CHANNEL_FREE, current, false, &curChannel);
		sounds->pop_back();
		channels->pop_back();
		ERRCHECK(result);

		// relase sounds
		// don't release yet, as soon as you release, the sound goes away
		//result = current->release();
		//ERRCHECK(result);
	}
	/*
	result = system->playSound(FMOD_CHANNEL_FREE, usersound, false, &channel1);
	ERRCHECK(result);
	*/
}

FMOD_RESULT F_CALLBACK pcmreadcallback(FMOD_SOUND *sound, void *data, unsigned int datalen)
{
    unsigned int  count;
    static float  t1 = 0, t2 = 0;        // time
    static float  v1 = 0, v2 = 0;        // velocity
    signed short *stereo16bitbuffer = (signed short *)data;

	/*
	for (count = 0; count<datalen; count++) {
		*stereo16bitbuffer++ = (signed short)(50 * sin((float)count));
	}
	*/

	
    for (count=0; count<datalen>>2; count++)        // >>2 = 16bit stereo (4 bytes per sample)
    {
        *stereo16bitbuffer++ = (signed short)(sin(t1) * 32767.0f);    // left channel
        *stereo16bitbuffer++ = (signed short)(sin(t2) * 32767.0f);    // right channel

        t1 += 0.01f   + v1;
        t2 += 0.0142f + v2;
        v1 += (float)(sin(t1) * 0.002f);
        v2 += (float)(sin(t2) * 0.002f);
    }

	

    return FMOD_OK; 
}


FMOD_RESULT F_CALLBACK pcmsetposcallback(FMOD_SOUND *sound, int subsound, unsigned int position, FMOD_TIMEUNIT postype)
{
    /*
        This is useful if the user calls Channel::setPosition and you want to seek your data accordingly.
    */
    return FMOD_OK;
}

//SoundManager * Sound_Manager = new SoundManager();



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
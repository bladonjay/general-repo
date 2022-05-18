/*
 * readNexFileC.c - A MEX function to read a NEX file (Nex Technologies),
 * and output the data in a MATLAB structure.
 *
 * NOTE: This script currently outputs raw tick marks and A/D values, use
 *       the accompanying wrapper script (readNexFile.m) to convert from
 *       tick marks to seconds and from A/D values to millivolts (mV).
 *       If you use readNexFileC.c to read the file, the format produced
 *       is the format expected by writeNexFileC.c. Similarly, if you use
 *       readNexFile.m to read the file, the format produced is the format
 *       expected by writeNexFile.m.
 *
 * USAGE:
 *   nex = readNexFileC('filename', verbose)
 *
 * INPUT:
 *   filename - Name of the NEX file to read.
 *   verbose  - Enable verbose output (optional, default 'false')
 *
 * OUTPUT:
 *   nex - A structure containing .nex file data
 *   nex.version - File version
 *   nex.comment - File comment
 *   nex.freq    - Frequency (to convert from ticks to seconds) (double)
 *   nex.tbeg    - Beginning of recording session (in ticks)    (int32)
 *   nex.tend    - End of recording session (in ticks)          (int32)
 *
 *   nex.neurons - Cell array of neuron structures
 *       neuron.name       - Name of neuron variable
 *       neuron.timestamps - Array of neuron timestamps (in ticks) (int32)
 *
 *   nex.events - Cell array of event structures
 *       event.name       - Name of event variable
 *       event.timestamps - Array of event timestamps (in ticks) (int32)
 *
 *   nex.intervals - Cell array of interval structures
 *       interval.name      - Name of interval variable
 *       interval.intStarts - Array of interval starts (in ticks) (int32)
 *       interval.intEnds   - Array of interval ends (in ticks) (int32)
 *
 *   nex.waves - Cell array of waveform structures
 *       wave.name        - Name of waveform variable
 *       wave.WFrequency  - A/D frequency for waveform data points (double)
 *       wave.ADtoMV      - Conversion from A/D values to mV (double)
 *       wave.NPointsWave - Number of data points in each wave (int32)
 *       wave.MVOffset    - Offset of A/D values (double)
 *                          (mv = raw*ADtoMV+MVOffset)
 *       wave.timestamps  - Array of waveform timestamps (in ticks) (int32)
 *       wave.waveforms   - Matrix of waveforms (in raw A/D values),
 *                          each column vector is a wave. (int16)
 *
 *   nex.popvectors - Cell array of population vector structures
 *       popvector.name    - Name of population vector variable
 *       popvector.weights - Array of population vector weights (double)
 *
 *   nex.contvars - Cell array of continuous variable structures
 *       contvar.name         - Name of continuous variable
 *       contvar.ADFrequency  - A/D frequency for data points (double)
 *       contvar.ADtoMV       - Conversion from A/D values to mV (double)
 *       contvar.MVOffset     - Offset of A/D values (double)
 *                              (mv = raw*ADtoMV+MVOffset)
 *
 *       continuous (A/D) data come in fragments. Each fragment has a
 *       timestamp and an index of the A/D data points in data array. The 
 *       timestamp corresponds to the time of recording of the first A/D
 *       value in this fragment.
 *
 *       contvar.timestamps     - Array of timestamps (in ticks) (int32)
 *       contvar.fragmentStarts - Array of start indexes (int32)
 *       contvar.data           - Array of data points (in raw A/D values) (int16)
 *
 *   nex.markers - Cell array of marker structures
 *       marker.name              - Name of marker variable
 *       marker.length            - Number of characters per marker value
 *       marker.timestamps        - Array of marker timestamps (in ticks)
 *       marker.values            - Array of marker value structures
 *           marker.value.name    - Name of marker value 
 *           marker.value.strings - Array of marker value strings
 *
 * Author: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
 * Last Modified: $Date: 2010-01-17 19:27:58 -0500 (Sun, 17 Jan 2010) $
 * Copyright (c) 2008-2010, Benjamin Kraus
 * $Id: readNexFileC.c 1107 2010-01-18 00:27:58Z bkraus $
 */

#include "mex.h"
#include "matrix.h"
#include "stdio.h"
#include "string.h"

/* The overall variable header should be 208 bytes. I added "MVOffset",
 * which is a double (8 bytes), but for some reason the structure size was
 * increasing by 16 bytes. Turns out it has to do with how the compiler
 * packs the data to increase efficiency. The "#pragma pack" command
 * changes the size used for packing, which resolves this problem.
 */
#pragma pack(push,4)

/* Nex file header structure */
struct NexFileHeader
{
	int    MagicNumber;    /* string NEX1 */
	int    Version;        /* 100 */
	char   Comment[256];
	double Frequency;      /* timestamped freq. - tics per second */
	int    Beg;            /* usually 0 */
	int    End;            /* = maximum timestamp + 1 */
	int    NumVars;        /* number of variables in the first batch */
	int    NextFileHeader; /* position of the next file header in the file
	                          not implemented yet */
	char   Padding[256];   /* future expansion */
};

/* Nex variable header structure */
struct NexVarHeader
{
	int    Type;         /* 0 - neuron, 1 event, 2- interval,
                          * 3 - waveform, 4 - population vector, 
                          * 5 - continuous variable, 6 - marker */
	int    Version;      /* 100 */
	char   Name[64];     /* variable name */
	int    DataOffset;   /* where the data array for this variable is located in the file */
	int    Count;        /* number of events, intervals, waveforms or weights */
	int    WireNumber;   /* neuron only, not used now */
	int    UnitNumber;   /* neuron only, not used now */
	int    Gain;         /* neuron only, not used now */
	int    Filter;       /* neuron only, not used now */
	double XPos;         /* neuron only, electrode position in (0,100) range, used in 3D */
	double YPos;         /* neuron only, electrode position in (0,100) range, used in 3D */
	double WFrequency;   /* waveform and continuous vars only, w/f sampling frequency */
	double ADtoMV;       /* waveform continuous vars only, coeff. to convert from A/D values to Millivolts. */
	int    NPointsWave;  /* waveform only, number of points in each wave */
	int    NMarkers;     /* how many values are associated with each marker */
	int    MarkerLength; /* how many characters are in each marker value */
	double MVOffset;     /* coeff to shift AD values in Millivolts: mv = raw*ADtoMV+MVOfffset */
	char   Padding[60];  /* header should be 208 bytes, add the necessary padding to make that true. */
};

#pragma pack(pop) /* Make sure to return packing to default setting. */

int readNexFile(mxArray *nexStruct, const char *fname, bool verbose);

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    int fnamelen, existans;
    char *fname;
    const char *fieldnames[12];
    
    bool verbose;

    mxArray *existin[2], *existout[1];

    if(nrhs < 1 || nrhs > 2)
        mexErrMsgTxt("Usage: nex = readNexFileC(filename, verbose)");
    else if(nlhs!=1)
        mexErrMsgTxt("Usage: nex = readNexFileC(filename, verbose)");
    else if(!mxIsChar(prhs[0]))
        mexErrMsgTxt("First argument must be a filename (string).");
    
    if(nrhs==1) verbose = false;
    else if(mxIsLogicalScalar(prhs[1]))
        verbose = mxIsLogicalScalarTrue(prhs[1]);
    else mexErrMsgTxt("Second argument must be logical true or false.");
    
    fnamelen = mxGetNumberOfElements(prhs[0])+1;
    fname = mxCalloc(fnamelen, sizeof(char));
    mxGetString(prhs[0],fname,fnamelen);
    
    fieldnames[0]  = "version";
    fieldnames[1]  = "comment";
    fieldnames[2]  = "freq";
    fieldnames[3]  = "tbeg";
    fieldnames[4]  = "tend";
    fieldnames[5]  = "neurons";
    fieldnames[6]  = "events";
    fieldnames[7]  = "intervals";
    fieldnames[8]  = "waves";
    fieldnames[9]  = "popvectors";
    fieldnames[10] = "contvars";
    fieldnames[11] = "markers";
    
    /* Base structure (nex) */
    plhs[0] = mxCreateStructMatrix(1,1,12,fieldnames);
    
    /* Check that the file is there to be read. To do this, call the matlab
     * function 'exist'. */
    existin[0] = mxCreateString(fname);
    existin[1] = mxCreateString("file");    
    mexCallMATLAB(1,existout,2,existin,"exist");
    
    existans = (int)mxGetScalar(existout[0]);
    switch (existans) {
        case 0: /* File/Directory not found. */
            mexPrintf("The file %s does not exist, aborting.\n",fname);
            break;
        case 2: /* File found, go ahead and read it. */
            readNexFile(plhs[0], fname, verbose);
            break;
        case 7: /* Directory found */
            mexPrintf("%s is a directory, aborting.\n",fname);
            break;
        default:
            mexPrintf("Unexpected output from 'exist', aborting.\n");
    }
}

int readNexFile(mxArray *nexStruct, const char *fname, bool verbose)
{
    FILE* fp;
    
	struct NexFileHeader fh;
    struct NexVarHeader *vh;
    char magic[] = "NEX1";
    
    const char *fieldnames[8][10];
    
    mxArray *neurons, *events, *intvals, *waves, *pops, *conts, *markers;
    mxArray *tmp_struct, *tmp, *tmp_cell, *values, *value_strs;
    
    int num_neurons = 0, num_events = 0, num_intvals = 0, num_waves = 0;
    int num_pops = 0, num_conts = 0, num_markers = 0, num_vars = 0;

    int cur_neuron = 0, cur_event = 0, cur_intval = 0, cur_wave = 0;
    int cur_pop = 0, cur_cont = 0, cur_marker = 0, cur_var = 0;
    
    int ntimes, nvals, markerlength, i, j, numread;
    
    char valname[64], **strbuf;
    
    fp = fopen(fname,"rb");
    if(fp == 0) {
        mexPrintf("Error opening file: %s\n", fname);
        return -1;
    }
   
    if(verbose) mexPrintf("Reading file header.\n");
    numread = fread(&fh, sizeof(struct NexFileHeader), 1, fp);
    if(numread != 1) {
        mexPrintf("Error reading file header: ");
        if(feof(fp)) mexPrintf("End of file reached early.\n");
        else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
        else mexPrintf("Unknown error.\n");
        fclose(fp); return -1;
    }
    
    if(fh.MagicNumber != *(int*)magic) {
        mexPrintf("File %s is not a valid NEX file.\n",fname);
        fclose(fp); return -1;
    }
    
    /* nex.version (int32) */
    tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    *(int *)mxGetData(tmp) = fh.Version;
    mxSetFieldByNumber(nexStruct,0, 0, tmp);
    
    /* nex.comment (string) */
    mxSetFieldByNumber(nexStruct,0, 1, mxCreateString(fh.Comment));
    
    /* nex.freq (double) */
    tmp = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    *mxGetPr(tmp) = fh.Frequency;
    mxSetFieldByNumber(nexStruct,0, 2, tmp);
    
    /* nex.tbeg (int32) */
    tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    *(int *)mxGetData(tmp) = fh.Beg;
    mxSetFieldByNumber(nexStruct,0, 3, tmp);
    
    /* nex.tend (int32) */
    tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    *(int *)mxGetData(tmp) = fh.End;
    mxSetFieldByNumber(nexStruct,0, 4, tmp);
    
    /* Read all the variable headers all at once. */
    vh = mxMalloc(fh.NumVars*sizeof(struct NexVarHeader));

    if(verbose) mexPrintf("Reading variable headers.\n");
    numread = fread(vh, sizeof(struct NexVarHeader), fh.NumVars, fp);
    if(numread != fh.NumVars) {
        mexPrintf("Error reading variables: ");
        if(feof(fp)) mexPrintf("End of file reached early.\n");
        else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
        else mexPrintf("Unknown error.\n");
        fclose(fp); return -1;
    }
    
    /* Tally the types of variable headers for use later. */
    num_vars = fh.NumVars;
    for (cur_var = 0; cur_var<num_vars; cur_var++) {
        switch (vh[cur_var].Type) {
            case 0: num_neurons++; break; /* Neuron              */
            case 1: num_events++;  break; /* Event               */
            case 2: num_intvals++; break; /* Interval            */
            case 3: num_waves++;   break; /* Wave                */
            case 4: num_pops++;    break; /* Population Vector   */
            case 5: num_conts++;   break; /* Continuous Variable */
            case 6: num_markers++; break; /* Marker              */            
        }
    }
    
    if(verbose) {
        mexPrintf("Total Variables: %d\n",num_vars);
        mexPrintf("Neurons: %d\n",num_neurons);
        mexPrintf("Events: %d\n",num_events);
        mexPrintf("Intervals: %d\n",num_intvals);
        mexPrintf("Waves: %d\n",num_waves);
        mexPrintf("Population Vectors: %d\n",num_pops);
        mexPrintf("Continuous Variables: %d\n",num_conts);
        mexPrintf("Markers: %d\n",num_markers);
    }
    
    /* Create the base cells for each variable type. */
    neurons = mxCreateCellMatrix(num_neurons,1);
    events  = mxCreateCellMatrix(num_events, 1);
    intvals = mxCreateCellMatrix(num_intvals,1);
    waves   = mxCreateCellMatrix(num_waves,  1);
    pops    = mxCreateCellMatrix(num_pops,   1);
    conts   = mxCreateCellMatrix(num_conts,  1);
    markers = mxCreateCellMatrix(num_markers,1);
    
    /* Assign the cells to the proper structure fields. */
    mxSetFieldByNumber(nexStruct,0, 5,neurons); /* Neurons               */
    mxSetFieldByNumber(nexStruct,0, 6,events);  /* Events                */
    mxSetFieldByNumber(nexStruct,0, 7,intvals); /* Intervals             */
    mxSetFieldByNumber(nexStruct,0, 8,waves);   /* Waveforms             */
    mxSetFieldByNumber(nexStruct,0, 9,pops);    /* Population vectors    */
    mxSetFieldByNumber(nexStruct,0,10,conts);   /* Continuous variables. */
    mxSetFieldByNumber(nexStruct,0,11,markers); /* Markers               */
    
    /* Prepare the names of the various structure fields for use later. */
    /* Neurons = Type 0 */
    fieldnames[0][0] = "name";
    fieldnames[0][1] = "version";
    fieldnames[0][2] = "wireNumber";
    fieldnames[0][3] = "unitNumber";
    fieldnames[0][4] = "gain";
    fieldnames[0][5] = "filter";
    fieldnames[0][6] = "xPos";
    fieldnames[0][7] = "yPos";
    fieldnames[0][8] = "timestamps";
    
    /* Events = Type 1 */
    fieldnames[1][0] = "name";
    fieldnames[1][1] = "version";
    fieldnames[1][2] = "timestamps";
    
    /* Intervals = Type 2 */
    fieldnames[2][0] = "name";
    fieldnames[2][1] = "version";
    fieldnames[2][2] = "intStarts";
    fieldnames[2][3] = "intEnds";
    
    /* Waves = Type 3 */
    fieldnames[3][0] = "name";
    fieldnames[3][1] = "version";
    fieldnames[3][2] = "wireNumber";
    fieldnames[3][3] = "unitNumber";
    fieldnames[3][4] = "WFrequency";
    fieldnames[3][5] = "ADtoMV";
    fieldnames[3][6] = "NPointsWave";
    fieldnames[3][7] = "MVOffset";
    fieldnames[3][8] = "timestamps";
    fieldnames[3][9] = "waveforms";
    
    /* Population Vectors = Type 4 */
    fieldnames[4][0] = "name";
    fieldnames[4][1] = "version";
    fieldnames[4][2] = "weights";
    
    /* Continuous Variables = Type 5 */
    fieldnames[5][0] = "name";
    fieldnames[5][1] = "version";
    fieldnames[5][2] = "ADFrequency";
    fieldnames[5][3] = "ADtoMV";
    fieldnames[5][4] = "MVOffset";
    fieldnames[5][5] = "timestamps";
    fieldnames[5][6] = "fragmentStarts";
    fieldnames[5][7] = "data";
    
    /* Markers = Type 6 */
    fieldnames[6][0] = "name";
    fieldnames[6][1] = "version";
    fieldnames[6][2] = "length";
    fieldnames[6][3] = "timestamps";
    fieldnames[6][4] = "values";
    
    /* Marker Values */
    fieldnames[7][0] = "name";
    fieldnames[7][1] = "strings";
    
    /* Start reading in the actual variable data. */
    if(verbose) mexPrintf("Reading variable data.\n");
    for (cur_var = 0; cur_var<num_vars; cur_var++) {
        switch (vh[cur_var].Type) {
            case 0: /* Neuron */
                /* Create the nex.neurons structure */
                tmp_struct = mxCreateStructMatrix(1,1,9,fieldnames[0]);
                mxSetCell(neurons, cur_neuron, tmp_struct);
                
                /* nex.neurons.name */
                mxSetFieldByNumber(tmp_struct,0,0,mxCreateString(vh[cur_var].Name));
                
                /* nex.neurons.version */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].Version;
                mxSetFieldByNumber(tmp_struct,0,1,tmp);
                
                /* nex.neurons.wireNumber */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].WireNumber;
                mxSetFieldByNumber(tmp_struct,0,2,tmp);
                
                /* nex.neurons.unitNumber */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].UnitNumber;
                mxSetFieldByNumber(tmp_struct,0,3,tmp);

                /* nex.neurons.gain */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].Gain;
                mxSetFieldByNumber(tmp_struct,0,4,tmp);
                
                /* nex.neurons.filter */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].Filter;
                mxSetFieldByNumber(tmp_struct,0,5,tmp);

                /* nex.neurons.xPos */
                tmp = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
                *mxGetPr(tmp) = vh[cur_var].XPos;
                mxSetFieldByNumber(tmp_struct,0,6,tmp);

                /* nex.neurons.yPos */
                tmp = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
                *mxGetPr(tmp) = vh[cur_var].YPos;
                mxSetFieldByNumber(tmp_struct,0,7,tmp);
                
                /* nex.neurons.timestamps */
                ntimes = vh[cur_var].Count;
                tmp = mxCreateNumericMatrix(ntimes,1,mxINT32_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,8,tmp);
                
                numread = fread((int *)mxGetData(tmp),sizeof(int),ntimes,fp);
                if(numread != ntimes) {
                    mexPrintf("Error reading neuron timestamps: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                cur_neuron++;
                break;
            case 1: /* Event */
                /* Create the nex.events structure */
                tmp_struct = mxCreateStructMatrix(1,1,3,fieldnames[1]);
                mxSetCell(events, cur_event, tmp_struct);
                
                /* nex.events.name */
                mxSetFieldByNumber(tmp_struct,0,0,mxCreateString(vh[cur_var].Name));
                
                /* nex.events.version */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].Version;
                mxSetFieldByNumber(tmp_struct,0,1,tmp);
                
                /* nex.events.timestamps */
                ntimes = vh[cur_var].Count;
                tmp = mxCreateNumericMatrix(ntimes,1,mxINT32_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,2,tmp);
                
                numread = fread((int *)mxGetData(tmp),sizeof(int),ntimes,fp);
                if(numread != ntimes) {
                    mexPrintf("Error reading event timestamps: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                cur_event++;
                break;
            case 2: /* Interval */
                /* Create the nex.intervals structure */
                tmp_struct = mxCreateStructMatrix(1,1,4,fieldnames[2]);
                mxSetCell(intvals, cur_intval, tmp_struct);
                
                /* nex.intervals.name */
                mxSetFieldByNumber(tmp_struct,0,0,mxCreateString(vh[cur_var].Name));
                
                /* nex.intervals.version */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].Version;
                mxSetFieldByNumber(tmp_struct,0,1,tmp);
                
                /* nex.intervals.intStarts */
                ntimes = vh[cur_var].Count;
                tmp = mxCreateNumericMatrix(ntimes,1,mxINT32_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,2,tmp);
                
                numread = fread((int *)mxGetData(tmp),sizeof(int),ntimes,fp);
                if(numread != ntimes) {
                    mexPrintf("Error reading interval start timestamps: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }                
                
                /* nex.intervals.intEnds */
                tmp = mxCreateNumericMatrix(ntimes,1,mxINT32_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,3,tmp);
                
                numread = fread((int *)mxGetData(tmp),sizeof(int),ntimes,fp);
                if(numread != ntimes) {
                    mexPrintf("Error reading interval end timestamps: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                cur_intval++;
                break;
            case 3: /* Wave */
                /* Create the nex.waves structure */
                tmp_struct = mxCreateStructMatrix(1,1,10,fieldnames[3]);
                mxSetCell(waves, cur_wave, tmp_struct);
                
                /* nex.waves.name */
                mxSetFieldByNumber(tmp_struct,0,0,mxCreateString(vh[cur_var].Name));
                
                /* nex.waves.version */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].Version;
                mxSetFieldByNumber(tmp_struct,0,1,tmp);
                
                /* nex.waves.wireNumber */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].WireNumber;
                mxSetFieldByNumber(tmp_struct,0,2,tmp);
                
                /* nex.waves.unitNumber */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].UnitNumber;
                mxSetFieldByNumber(tmp_struct,0,3,tmp);

                /* nex.waves.WFrequency */
                tmp = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
                *mxGetPr(tmp) = vh[cur_var].WFrequency;
                mxSetFieldByNumber(tmp_struct,0,4,tmp);
                
                /* nex.waves.ADtoMV */
                tmp = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
                *mxGetPr(tmp) = vh[cur_var].ADtoMV;
                mxSetFieldByNumber(tmp_struct,0,5,tmp);
                
                /* nex.waves.NPointsWave */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].NPointsWave;
                mxSetFieldByNumber(tmp_struct,0,6,tmp);
                
                /* nex.waves.MVOffset */
                tmp = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
                *mxGetPr(tmp) = vh[cur_var].MVOffset;
                mxSetFieldByNumber(tmp_struct,0,7,tmp);
                
                /* nex.waves.timestamps */
                ntimes = vh[cur_var].Count;
                tmp = mxCreateNumericMatrix(ntimes,1,mxINT32_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,8,tmp);
                
                numread = fread((int *)mxGetData(tmp),sizeof(int),ntimes,fp);
                if(numread != ntimes) {
                    mexPrintf("Error reading wave timestamps: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                /* nex.waves.waveforms */
                nvals = vh[cur_var].NPointsWave;
                tmp = mxCreateNumericMatrix(nvals,ntimes,mxINT16_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,9,tmp);
                
                numread = fread((short *)mxGetData(tmp),sizeof(short),ntimes*nvals,fp);                
                if(numread != ntimes*nvals) {
                    mexPrintf("Error reading wave data: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                cur_wave++;
                break;
            case 4: /* Population Vector */
                /* Create the nex.popvectors structure */
                tmp_struct = mxCreateStructMatrix(1,1,3,fieldnames[4]);
                mxSetCell(pops, cur_pop, tmp_struct);
                
                /* nex.popvectors.name */
                mxSetFieldByNumber(tmp_struct,0,0,mxCreateString(vh[cur_var].Name));
                
                /* nex.popvectors.version */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].Version;
                mxSetFieldByNumber(tmp_struct,0,1,tmp);
                
                /* nex.popvectors.weights */
                ntimes = vh[cur_var].Count;
                tmp = mxCreateNumericMatrix(ntimes,1,mxDOUBLE_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,2,tmp);
                
                numread = fread(mxGetPr(tmp),sizeof(double),ntimes,fp);
                if(numread != ntimes) {
                    mexPrintf("Error reading population vector weights: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                cur_pop++;
                break;
            case 5: /* Continuous Variable */
                /* Create the nex.contvars structure */
                tmp_struct = mxCreateStructMatrix(1,1,8,fieldnames[5]);
                mxSetCell(conts, cur_cont, tmp_struct);
                
                /* nex.contvars.name */
                mxSetFieldByNumber(tmp_struct,0,0,mxCreateString(vh[cur_var].Name));
                
                /* nex.contvars.version */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].Version;
                mxSetFieldByNumber(tmp_struct,0,1,tmp);
                
                /* nex.contvars.ADFrequency */
                tmp = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
                *mxGetPr(tmp) = vh[cur_var].WFrequency;
                mxSetFieldByNumber(tmp_struct,0,2,tmp);
                
                /* nex.contvars.ADtoMV */
                tmp = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
                *mxGetPr(tmp) = vh[cur_var].ADtoMV;
                mxSetFieldByNumber(tmp_struct,0,3,tmp);
                
                /* nex.contvars.MVOffset */
                tmp = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
                *mxGetPr(tmp) = vh[cur_var].MVOffset;
                mxSetFieldByNumber(tmp_struct,0,4,tmp);
                
                /* nex.contvars.timestamps */
                ntimes = vh[cur_var].Count;
                tmp = mxCreateNumericMatrix(ntimes,1,mxINT32_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,5,tmp);
                
                numread = fread((int *)mxGetData(tmp),sizeof(int),ntimes,fp);
                if(numread != ntimes) {
                    mexPrintf("Error reading continous variable timestamps: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                /* nex.contvars.fragmentStarts */
                tmp = mxCreateNumericMatrix(ntimes,1,mxINT32_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,6,tmp);
                
                numread = fread((int *)mxGetData(tmp),sizeof(int),ntimes,fp);
                if(numread != ntimes) {
                    mexPrintf("Error reading continuous variable fragment indexes: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                /* nex.contvars.data */
                nvals = vh[cur_var].NPointsWave;
                tmp = mxCreateNumericMatrix(nvals,1,mxINT16_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,7,tmp);
                
                numread = fread((short *)mxGetData(tmp),sizeof(short),nvals,fp);
                if(numread != nvals) {
                    mexPrintf("Error reading continuous variable data: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                cur_cont++;
                break;
            case 6: /* Marker */
                /* Create the nex.markers structure */
                tmp_struct = mxCreateStructMatrix(1,1,5,fieldnames[6]);
                mxSetCell(markers, cur_marker, tmp_struct);
                
                /* nex.markers.name */
                mxSetFieldByNumber(tmp_struct,0,0,mxCreateString(vh[cur_var].Name));
                
                /* nex.markers.version */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].Version;
                mxSetFieldByNumber(tmp_struct,0,1,tmp);
                
                /* nex.markers.length */
                tmp = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                *(int *)mxGetData(tmp) = vh[cur_var].MarkerLength;
                mxSetFieldByNumber(tmp_struct,0,2,tmp);
                
                /* nex.markers.timestamps */
                ntimes = vh[cur_var].Count;
                tmp = mxCreateNumericMatrix(ntimes,1,mxINT32_CLASS,mxREAL);
                mxSetFieldByNumber(tmp_struct,0,3,tmp);
                
                numread = fread((int *)mxGetData(tmp),sizeof(int),ntimes,fp);
                if(numread != ntimes) {
                    mexPrintf("Error reading marker timestamps: ");
                    if(feof(fp)) mexPrintf("End of file reached early.\n");
                    else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                    else mexPrintf("Unknown error.\n");
                    fclose(fp); return -1;
                }
                
                /* nex.markers.values */
                nvals = vh[cur_var].NMarkers;
                tmp_cell = mxCreateCellMatrix(nvals,1);
                mxSetFieldByNumber(tmp_struct,0,4,tmp_cell);
                
                markerlength = vh[cur_var].MarkerLength;
                strbuf = mxMalloc(ntimes*sizeof(char *));
                
                for (i = 0; i<nvals; i++) {
                    /* nex.markers.values */
                    values = mxCreateStructMatrix(1,1,2,fieldnames[7]);
                    mxSetCell(tmp_cell, i, values);
                    
                    /* nex.markers.values.name */
                    numread = fread(valname,sizeof(char),64,fp);
                    if(numread != 64) {
                        mexPrintf("Error reading marker value names: ");
                        if(feof(fp)) mexPrintf("End of file reached early.\n");
                        else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                        else mexPrintf("Unknown error.\n");
                        fclose(fp); return -1;
                    }
                    
                    mxSetFieldByNumber(values,0,0,mxCreateString(valname));
                    
                    /* nex.markers.values.strings */
                    for (j = 0; j<ntimes; j++) {
                        strbuf[j] = mxMalloc(markerlength*sizeof(char));
                        numread = fread(strbuf[j],sizeof(char),markerlength,fp);
                        if(numread != markerlength) {
                            mexPrintf("Error reading marker value strings: ");
                            if(feof(fp)) mexPrintf("End of file reached early.\n");
                            else if(ferror(fp)) mexPrintf("Error code %d.\n",ferror(fp));
                            else mexPrintf("Unknown error.\n");
                            fclose(fp); return -1;
                        }
                    }
                    if(ferror(fp)) {
                        mexPrintf("Error reading marker strings: ");
                        if(feof(fp)) mexPrintf("End of file reached early.\n");
                        else mexPrintf("Error code %d.\n",ferror(fp));
                        fclose(fp); return -1;
                    }
                    
                    value_strs = mxCreateCharMatrixFromStrings(ntimes, (const char **)strbuf);
                    
                    mxSetFieldByNumber(values,0,1,value_strs);
                    
                    for (j = 0; j<ntimes; j++) mxFree(strbuf[j]);
                }
                mxFree(strbuf);
                
                cur_marker++;
                break;
        }
    }
        
    fclose(fp);
    return 0;
}


#ifndef FASTQ_COMMONS_H
#define FASTQ_COMMONS_H

#define VALID_READS_FILE_SUFFIX 	".valid"
#define INVALID_READS_FILE_SUFFIX 	".invalid"

#define	ALL_READS			1
#define	ONLY_VALID_READS		2
#define	ONLY_INVALID_READS		3

#define MIN_QUALITY_VALUE		10	//Normalized
#define MAX_QUALITY_VALUE		70	//Normalized

/* **********************************************
 *  	Global variables and structures		*
 * *********************************************/

extern double chaos_game_time;
extern struct timeval t1_chaos_game, t2_chaos_game;

#endif /* FASTQ_COMMONS_H */
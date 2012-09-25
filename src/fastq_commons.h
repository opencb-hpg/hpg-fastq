/*
 * Copyright (c) 2012 Victor Requena (BULL)
 * Copyright (c) 2012 Ignacio Medina (CGI-CIPF)
 *
 * This file is part of hpg-fastq.
 *
 * hpg-fastq is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hpg-fastq is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-fastq. If not, see <http://www.gnu.org/licenses/>.
 */


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
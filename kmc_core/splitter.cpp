/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.1.1
Date   : 2019-05-19
*/

#include "splitter.h"
#include <iostream>

//************************************************************************************************************
// CSplitter class - splits kmers into bins according to their signatures
//************************************************************************************************************

//uint32 CSplitter::MAX_LINE_SIZE = 1 << 14;
uint32 CSplitter::MAX_LINE_SIZE = 1 << 16;

//----------------------------------------------------------------------------------
// Assigns queues
CSplitter::CSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	//mm = Queues.mm;
	file_type = Params.file_type;
	both_strands = Params.both_strands;

	bin_part_queue = Queues.bpq.get();
	pmm_reads = Queues.pmm_reads.get();
	kmer_len = Params.kmer_len;
	signature_len = Params.signature_len;

	mem_part_pmm_bins = Params.mem_part_pmm_bins;

	mem_part_pmm_reads = Params.mem_part_pmm_reads; 

	s_mapper = Queues.s_mapper.get();

	part = nullptr;

	// Prepare encoding of symbols
	for (int i = 0; i < 256; ++i)
		codes[i] = -1;
	codes['A'] = codes['a'] = 0;
	codes['C'] = codes['c'] = 1;
	codes['G'] = codes['g'] = 2;
	codes['T'] = codes['t'] = 3;

	n_reads = 0;

	homopolymer_compressed = Params.homopolymer_compressed;

	ntHashEstimator = Queues.ntHashEstimator.get();
}

void CSplitter::InitBins(CKMCParams &Params, CKMCQueues &Queues)
{
	n_bins = Params.n_bins;
	uint32 buffer_size = Params.bin_part_size;
	// Create objects for all bin
	bins.resize(n_bins);
	for (uint32 i = 0; i < n_bins; ++i)
	{
		bins[i] = std::make_unique<CKmerBinCollector>(Queues, Params, buffer_size, i);
	}
}

//----------------------------------------------------------------------------------
// Parse long read, header_merker is '@' or '>'
bool CSplitter::GetSeqLongRead(char *seq, uint32 &seq_size, uchar header_marker)
{	
	uint32 pos = 0;
	//long read may or may not contain header
	if (part_pos == 0 && part[0] == header_marker)
	{
		++n_reads;
		for (; part[part_pos] != '\n' && part[part_pos] != '\r'; ++part_pos)
			;
	}
	while (pos < mem_part_pmm_reads && part_pos < part_size)
		seq[pos++] = codes[part[part_pos++]];
	seq_size = pos;
	if (part_pos < part_size)
		part_pos -= kmer_len - 1;
	return true;
}

//----------------------------------------------------------------------------------
// Return a single record from FASTA/FASTQ data
bool CSplitter::GetSeq(char *seq, uint32 &seq_size, ReadType read_type)
{
	if (part_pos >= part_size)
		return false;

	uchar c = 0;
	uint32 pos = 0;

	if (file_type == InputType::FASTA || file_type == InputType::KMC)
	{		
		if (read_type == ReadType::long_read)
			return GetSeqLongRead(seq, seq_size, '>');
		if (curr_read_len == 0)
		{
			// Title
			c = part[part_pos++];
			if (c != '>')
				return false;
			++n_reads;

			for (; part_pos < part_size;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
			}
			if (part_pos >= part_size)
				return false;

			c = part[part_pos++];
			if (c >= 32 || c == part[part_pos - 2]) //read may be empty
				part_pos--;
			else if (part_pos >= part_size)
				return false;

			// Sequence
			for (; part_pos < part_size && pos < mem_part_pmm_reads;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
				seq[pos++] = codes[c];
			}

			seq_size = pos;

			if (part_pos >= part_size)
				return true;

			curr_read_len = pos;

			if (pos >= mem_part_pmm_reads) // read is too long to fit into out buff, it will be splitted into multiple buffers
			{
				part_pos -= kmer_len - 1;
				return true;
			}
		}
		else // we are inside read
		{
			// Sequence
			for (; part_pos < part_size && pos < mem_part_pmm_reads;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
				seq[pos++] = codes[c];
			}

			seq_size = pos;

			if (part_pos >= part_size)
				return true;

			curr_read_len += pos - kmer_len + 1;

			if (pos >= mem_part_pmm_reads) // read is too long to fit into out buff, it will be splitted into multiple buffers
			{
				part_pos -= kmer_len - 1;
				return true;
			}
		}
		
		curr_read_len = 0;

		//end of last record 
		if (part_pos >= part_size)
			return true;

		if (part[part_pos++] >= 32)
			part_pos--;
		else if (part_pos >= part_size)
			return true;
	}
	else if (file_type == InputType::FASTQ)
	{
		if (read_type == ReadType::long_read)		
			return GetSeqLongRead(seq, seq_size, '@');	
	
		if (curr_read_len == 0)
		{
			// Title
			c = part[part_pos++];
			if (c != '@')
				return false;
			++n_reads;

			for (; part_pos < part_size;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
			}
			if (part_pos >= part_size)
				return false;

			c = part[part_pos++];
			if (c >= 32 || c == part[part_pos - 2]) //read may be empty
				part_pos--;
			else if (part_pos >= part_size)
				return false;

			// Sequence
			for (; part_pos < part_size && pos < mem_part_pmm_reads;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
				seq[pos++] = codes[c];
			}
			if (part_pos >= part_size)
				return false;

			seq_size = pos;
			curr_read_len = pos;

			if (pos >= mem_part_pmm_reads) // read is too long to fit into out buff, it will be splitted into multiple buffers
			{
				part_pos -= kmer_len - 1;
				return true;
			}
		}
		else // we are inside read
		{
			// Sequence
			for (; part_pos < part_size && pos < mem_part_pmm_reads;)
			{
				c = part[part_pos++];
				if (c < 32)					// newliners
					break;
				seq[pos++] = codes[c];
			}
			if (part_pos >= part_size)
				return false;

			seq_size = pos;
			curr_read_len += pos - kmer_len + 1;
			if (pos >= mem_part_pmm_reads) // read is too long to fit into out buff, it will be splitted into multiple buffers
			{
				part_pos -= kmer_len - 1;
				return true;
			}
		}

		c = part[part_pos++];
		if (c >= 32)
			part_pos--;
		else if (part_pos >= part_size)
			return false;

		// Plus
		c = part[part_pos++];
		if (part_pos >= part_size)
			return false;
		if (c != '+')
			return false;
		for (; part_pos < part_size;)
		{
			c = part[part_pos++];
			if (c < 32)					// newliners
				break;
		}
		if (part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if (c >= 32 || c == part[part_pos - 2]) //qual may be empty
			part_pos--;
		else if (part_pos >= part_size)
			return false;

		// Quality
		part_pos += curr_read_len; //skip quality
		
		curr_read_len = 0;

		if (part_pos >= part_size)
			return false;
		c = part[part_pos++];

		//end of last record 
		if (part_pos >= part_size)
			return true;

		//may be additional EOL character 
		if (part[part_pos++] >= 32)
			part_pos--;
		else if (part_pos >= part_size)
			return true;
		
	}
	else if (file_type == InputType::MULTILINE_FASTA)
	{
		if (part[part_pos] == '>')//need to ommit header
		{
			++n_reads;
			for (; part_pos < part_size && part[part_pos] != '\n' && part[part_pos] != '\r'; ++part_pos);//find EOF
			++part_pos;
			if (part[part_pos] == '\n' || part[part_pos] == '\r')
				++part_pos;
		}
		for (; part_pos < part_size && pos < mem_part_pmm_reads && part[part_pos] != '>';)
		{
			seq[pos++] = codes[part[part_pos++]];
		}
		seq_size = pos;
		if (part_pos < part_size && part[part_pos] != '>')//need to copy last k-1 kmers 
		{
			part_pos -= kmer_len - 1;
		}
		return true;

	}
	else if (file_type == InputType::BAM)
	{
		while (true)
		{
			if (part_pos >= part_size)
				return false;

			int32_t block_size;
			read_int32_t(block_size, part, part_pos);

			uint64_t start_pos = part_pos;

			part_pos += 8;

			uint32_t bin_mq_nl;
			read_uint32_t(bin_mq_nl, part, part_pos);

			uint32_t l_read_name = (bin_mq_nl & ((1 << 8) - 1));
			uint32_t flag_nc;
			read_uint32_t(flag_nc, part, part_pos);
			uint32_t n_cigar_op = flag_nc & ((1ul << 16) - 1);
			int32_t l_seq;
			read_int32_t(l_seq, part, part_pos);

			part_pos += 12;

			uint32_t flags = flag_nc >> 16;

			bool exclude_read = ((flags >> 8) & 1) || ((flags >> 11) & 1); //TODO: I think that is the way samtools filter out some reads (secondary and supplementary)

			part_pos += l_read_name; // skip read name

			part_pos += 4 * n_cigar_op;
			if (!exclude_read)
			{
				bool is_rev_comp = (flags >> 4) & 1;

				if (!both_strands && is_rev_comp) //if read is reversed and kmc was run to count all (not only canonical) kmers read must be transformed back
				{
					//static const char rev_maping[] = "=TGMCRSVAWYHKDBN";
					static const char rev_maping[] = { -1, 3, 2, -1, 1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1 };// "=TGMCRSVAWYHKDBN";
					uint32 n_bytes = l_seq / 2;
					uint64_t pos_after = pos + l_seq;
					pos = pos_after;
					for (uint32_t ii = 0; ii < n_bytes; ++ii)
					{
						unsigned char byte = part[part_pos++];
						seq[--pos_after] = rev_maping[byte >> 4];
						seq[--pos_after] = rev_maping[byte & 15];
					}

					if (l_seq & 1) //odd
					{
						unsigned char byte = part[part_pos++];
						seq[--pos_after] = rev_maping[byte >> 4];
					}
				}
				else
				{
					static const char maping[] = { -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1 };//"=ACMGRSVTWYHKDBN";
					uint32 n_bytes = l_seq / 2;
					for (uint32_t ii = 0; ii < n_bytes; ++ii)
					{
						unsigned char byte = part[part_pos++];
						seq[pos++] = maping[byte >> 4];
						seq[pos++] = maping[byte & 15];
					}

					if (l_seq & 1) //odd
					{
						unsigned char byte = part[part_pos++];
						seq[pos++] = maping[byte >> 4];
					}
				}
				seq_size = pos;
			}
			else
			{
				part_pos += (l_seq + 1) / 2;
			}

			//move to next record
			uint64_t readed = part_pos - start_pos;
			uint64_t remaining = block_size - readed;
			part_pos += remaining;

			if (!exclude_read) //if readed successfuly return		
			{
				++n_reads;
				return true;
			}
		}

	}
	return (c == '\n' || c == '\r');
}
//----------------------------------------------------------------------------------
void CSplitter::HomopolymerCompressSeq(char* seq, uint32 &seq_size)
{
	if (seq_size <= 1)
		return;

	uint32 read_pos = 1;
	uint32 write_pos = 0;
	for (; read_pos < seq_size; ++read_pos)
		if (seq[read_pos] != seq[write_pos])
			seq[++write_pos] = seq[read_pos];
	seq_size = write_pos + 1;
}

//----------------------------------------------------------------------------------
// Calculate statistics of m-mers
void CSplitter::CalcStats(uchar* _part, uint64 _part_size, ReadType read_type, uint32* _stats)
{

	part = _part;
	part_size = _part_size;
	part_pos = 0;

	char *seq;
	uint32 seq_size;
	pmm_reads->reserve(seq);
	CMmer current_signature(signature_len);
	uint32 i;

	// ADDED VARIABLES
	uint32_t w_len = kmer_len; // Window length. Using w=k for now
	uint32_t canonical_flag = 1; // 1 for canonical mode and 0 for forward strand only
	uint32_t min_flag = 0; // checks if a minimizer has already been computed in an iteration
	uint64_t kmer_int = 0; // Integer representation of a kmer
	uint64_t rcm_kmer_int = 0; // Integer representation for the reverse complement of a kmer
	uint64_t can_int = 0; // Stores smaller of hash(kmer_int) and hash(rcm_kmer_int)
	uint64_t kmer_strand = 0; // 1 if the kmer was from the reverse strand and 0 otherwise
	uint64_t mask1 = (1ULL<<2 * kmer_len) - 1; // Mask to keep the kmer_int values in range
  	uint64_t shift1 = 2 * (kmer_len-1);
	uint64_t kmer_hash = UINT64_MAX; // Stores the hash value of the kmer formed
	uint64_t min_hash = UINT64_MAX; // Stores the hash value of the last minimizer
	uint64_t min_pos = 0; // Stores index of the last minimizer in the buffer
	uint64_t min_loc = 0; // Stores the index of the last minimizer in the sequence
	uint64_t min_strand = 0; // 1 if minimizer was from reverse strand and 0 otherwise
	uint64_t buf[256]; // buffer to store w kmer_hashes at a time
	uint64_t buf_pos = 0; // index for the buffer
	uint64_t kmer_strand_buf[256]; // buffer to store w kmer_strand values at a time
	char *rev; // stores the reverse complement of a kmer
	pmm_reads->reserve(rev);

	while (GetSeq(seq, seq_size, read_type))
	{

		if (homopolymer_compressed)
			HomopolymerCompressSeq(seq, seq_size);

		i = 0;
		bool contains_N = false;
		while (i + kmer_len - 1 < seq_size)
		{		
			for (; i < seq_size; ++i)
			{
				if(seq[i] < 0)
				{
					contains_N = true;
					break;
				}
				else
				{
					kmer_int = (kmer_int << 2 | seq[i]) & mask1;
					kmer_strand = 0;
					can_int = kmer_int;
					if (canonical_flag) 
					{
						rcm_kmer_int = (rcm_kmer_int >> 2) | (3ULL^seq[i]) << shift1;
						if(rcm_kmer_int < kmer_int){ // rcm hash is smaller so use that 
							can_int = rcm_kmer_int;
							kmer_strand = 1;
						}
					}
				}
				
				if(i>=kmer_len-1) // atleast one full kmer formed
				{
					kmer_hash = hash64(can_int, mask1) << 8 | kmer_len;
					buf[buf_pos] = kmer_hash; // push the kmer hash to the buffer
					kmer_strand_buf[buf_pos] = kmer_strand;
				}

				if(kmer_hash <= min_hash && i<w_len+kmer_len-1 && i>=kmer_len-1) // get right most minimizer in first window
				{
					min_hash = kmer_hash;
					min_pos = buf_pos;
					min_loc = i;
					min_flag = 1;
					min_strand = kmer_strand;
				}

				if(kmer_hash < min_hash && i>=w_len+kmer_len-1) // smaller kmer found
				{
					if(min_hash != UINT64_MAX)
					{
						// add the current minimizer
						
						// Minimizer came from forward strand so form the signature as is
						if(min_strand == 0){
							current_signature.insert(seq + min_loc - kmer_len + 1);
							_stats[current_signature.get()] += 1;
						}

						// Minimizer came from reverse strand. SO rev comp and then form the signature
						// Ensures that both the kmer and its reverse complement have the same signature so that they are added to the same bins
						else{
							
							for(uint32_t rev_i = 0; rev_i < signature_len; rev_i++){
								rev[rev_i] = 3 - (int)seq[min_loc - rev_i];
							}
							current_signature.insert(rev);
							_stats[current_signature.get()] += 1;
						}
					}

					min_hash = kmer_hash;
					min_pos = buf_pos;
					min_loc = i;
					min_flag = 1;
					min_strand = kmer_strand;
				}

				if(buf_pos == min_pos && min_flag == 0 && i>=kmer_len-1) // old minimizer moved out of window
				{

					if(i >= w_len + kmer_len - 1 && min_hash != UINT64_MAX)
					{
						// add the current minimizer
						if(min_strand == 0){
							current_signature.insert(seq + min_loc - kmer_len + 1);
							_stats[current_signature.get()] += 1;
						}
						else{
							
							for(uint32_t rev_i = 0; rev_i < signature_len; rev_i++){
								rev[rev_i] = 3 - (int)seq[min_loc - rev_i];
							}
							current_signature.insert(rev);
							_stats[current_signature.get()] += 1;
						}
					}

					min_hash = UINT64_MAX;

					uint32_t j;
					uint32_t j_counter = 0;

					// Loop over the past window to find the rightmost minimizer

					for (j = buf_pos + 1; j < w_len; ++j)
					{ 
						if (buf[j] <= min_hash) //  >= is important s.t. min is always the closest k-mer
						{
							min_hash = buf[j];
							min_pos = j;
							min_loc = i - w_len + 1 + j_counter; // min_loc is j_counter times ahead of the start of the buffer
							min_strand = kmer_strand_buf[j];
						} 
						
						j_counter++;
					}

					for (j = 0; j <= buf_pos; ++j)
					{
						if (buf[j] <= min_hash)
						{
							min_hash = buf[j];
							min_pos = j;
							min_loc = i - w_len + 1 + j_counter;
							min_strand = kmer_strand_buf[j];
						}

						j_counter++;
					}

				}

				min_flag = 0; // reset min_flag to 0 for next run
				if(i>=kmer_len-1)
				{
					if (++buf_pos == w_len) buf_pos = 0; // if buffer is full overwrite the oldest kmer
				}
			}

			// If a read contains an N character ignore the read from there on
			if (contains_N)
			{
				break;
			}
		}

		if(!contains_N)
		{
			// adding the last minimizer in the read

			if(min_strand == 0){
				current_signature.insert(seq + min_loc - kmer_len + 1);
				_stats[current_signature.get()] += 1;
			}
			else{
				for(uint32_t rev_i = 0; rev_i < signature_len; rev_i++){
					rev[rev_i] = 3 - (int)seq[min_loc - rev_i];
				}
				
				current_signature.insert(rev);
				_stats[current_signature.get()] += 1;
			}

		}
	}

	pmm_reads->free(seq);
	pmm_reads->free(rev);

}

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part, but only for k-mer occurence estimation
bool CSplitter::ProcessReadsOnlyEstimate(uchar* _part, uint64 _part_size, ReadType read_type)
{
	part = _part;
	part_size = _part_size;
	part_pos = 0;

	char* seq;
	uint32 seq_size;
	pmm_reads->reserve(seq);

	while (GetSeq(seq, seq_size, read_type))
		ntHashEstimator->Process(seq, seq_size);

	pmm_reads->free(seq);

	return true;
}

// -------------------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part
bool CSplitter::ProcessReads(uchar *_part, uint64 _part_size, ReadType read_type)
{
	part = _part;
	part_size = _part_size;
	part_pos = 0;

	char *seq;
	uint32 seq_size;
	pmm_reads->reserve(seq);
	CMmer current_signature(signature_len);
	uint32 bin_no;
	uint32 i;

	// ADDED VARIABLES
	uint32_t w_len = kmer_len; // Window length. Using w=k for now
	uint32_t canonical_flag = 1; // 1 for canonical mode and 0 for forward strand only
	uint32_t min_flag = 0; // checks if a minimizer has already been computed in an iteration
	uint64_t kmer_int = 0; // Integer representation of a kmer
	uint64_t rcm_kmer_int = 0; // Integer representation for the reverse complement of a kmer
	uint64_t can_int = 0; // Stores smaller of hash(kmer_int) and hash(rcm_kmer_int)
	uint64_t kmer_strand = 0; // 1 if the kmer was from the reverse strand and 0 otherwise
	uint64_t mask1 = (1ULL<<2 * kmer_len) - 1; // Mask to keep the kmer_int values in range
  	uint64_t shift1 = 2 * (kmer_len-1);
	uint64_t kmer_hash = UINT64_MAX; // Stores the hash value of the kmer formed
	uint64_t min_hash = UINT64_MAX; // Stores the hash value of the last minimizer
	uint64_t min_pos = 0; // Stores index of the last minimizer in the buffer
	uint64_t min_loc = 0; // Stores the index of the last minimizer in the sequence
	uint64_t min_strand = 0; // 1 if minimizer was from reverse strand and 0 otherwise
	uint64_t buf[256]; // buffer to store w kmer_hashes at a time
	uint64_t buf_pos = 0; // index for the buffer
	uint64_t kmer_strand_buf[256]; // buffer to store w kmer_strand values at a time
	char *rev; // stores the reverse complement of a kmer
	pmm_reads->reserve(rev);
	

	while (GetSeq(seq, seq_size, read_type))
	{
		if (ntHashEstimator)
			ntHashEstimator->Process(seq, seq_size);

		if (homopolymer_compressed)
			HomopolymerCompressSeq(seq, seq_size);
		
		i = 0;
		bool contains_N = false;
		while (i + kmer_len - 1 < seq_size)
		{
			for (; i < seq_size; ++i)
			{
				if(seq[i] < 0)
				{
					++i;
					contains_N = true;
					break;
				}
				else
				{
					kmer_int = (kmer_int << 2 | seq[i]) & mask1;
					kmer_strand = 0;
					can_int = kmer_int;
					if (canonical_flag) 
					{
						rcm_kmer_int = (rcm_kmer_int >> 2) | (3ULL^seq[i]) << shift1;
						if(rcm_kmer_int < kmer_int){
							can_int = rcm_kmer_int;
							kmer_strand = 1;
						}
					}
				}

				if(i>=kmer_len-1) // atleast one full kmer formed
				{
					kmer_hash = hash64(can_int, mask1) << 8 | kmer_len;
					buf[buf_pos] = kmer_hash;
					kmer_strand_buf[buf_pos] = kmer_strand;
				}
				
				if(kmer_hash <= min_hash && i<w_len+kmer_len-1 && i>=kmer_len-1) // get right most in first window
				{
					min_hash = kmer_hash;
					min_pos = buf_pos;
					min_loc = i;
					min_flag = 1;
					min_strand = kmer_strand;
				}

				if(kmer_hash < min_hash && i>=w_len+kmer_len-1) // smaller kmer found
				{
					if(min_hash != UINT64_MAX)
					{
						// add the current minimizer
						
						// Minimizer came from forward strand so form the signature as is
						if(min_strand == 0){
							current_signature.insert(seq + min_loc - kmer_len + 1);
							bin_no = s_mapper->get_bin_id(current_signature.get()); // compute bin for the minimizer using the signature
							bins[bin_no]->PutExtendedKmer(seq + min_loc - kmer_len + 1, kmer_len); // push the minimizer to the allocated bin
						}

						// Minimizer came from reverse strand. SO rev comp and then form the signature
						// Ensures that both the kmer and its reverse complement have the same signature so that they are added to the same bins
						else{
							
							for(uint32_t rev_i = 0; rev_i < signature_len; rev_i++){
								rev[rev_i] = 3 - (int)seq[min_loc - rev_i];
							}
							current_signature.insert(rev);
							bin_no = s_mapper->get_bin_id(current_signature.get()); 
							bins[bin_no]->PutExtendedKmer(seq + min_loc - kmer_len + 1, kmer_len);
						}
					}

					min_hash = kmer_hash;
					min_pos = buf_pos;
					min_loc = i;
					min_flag = 1;
					min_strand = kmer_strand;
				}

				if(buf_pos == min_pos && min_flag == 0 && i>=kmer_len-1) // old minimizer moved out of window
				{
					if(i >= w_len + kmer_len - 1 && min_hash != UINT64_MAX)
					{
						// add the current minimizer
						if(min_strand == 0){
							current_signature.insert(seq + min_loc - kmer_len + 1);
							bin_no = s_mapper->get_bin_id(current_signature.get());
							bins[bin_no]->PutExtendedKmer(seq + min_loc - kmer_len + 1, kmer_len);
						}
						else{
							for(uint32_t rev_i = 0; rev_i < signature_len; rev_i++){
								rev[rev_i] = 3 - (int)seq[min_loc - rev_i];
							}

							current_signature.insert(rev);
							bin_no = s_mapper->get_bin_id(current_signature.get());
							bins[bin_no]->PutExtendedKmer(seq + min_loc - kmer_len + 1, kmer_len);
						}
					}

					min_hash = UINT64_MAX;

					uint32_t j;
					uint32_t j_counter = 0;

					// Loop over the past window to find the rightmost minimizer
					
					for (j = buf_pos + 1; j < w_len; ++j)
					{ 
						if (buf[j] <= min_hash) //  >= is important s.t. min is always the closest k-mer
						{
							min_hash = buf[j];
							min_strand = kmer_strand_buf[j];
							min_pos = j;
							min_loc = i - w_len + 1 + j_counter;
						} 
						
						j_counter++;
					}
					for (j = 0; j <= buf_pos; ++j)
					{
						if (buf[j] <= min_hash)
						{
							min_hash = buf[j];
							min_strand = kmer_strand_buf[j];
							min_pos = j;
							min_loc = i - w_len + 1 + j_counter;
						}

						j_counter++;
					}
				}

				min_flag = 0;
				if(i>=kmer_len-1)
				{
					if (++buf_pos == w_len) buf_pos = 0; // if buffer is full overwrite the oldest kmer
				}
				
			}
			if (contains_N)
			{
				break;
			}
		}

		if(!contains_N){

			// adding the last minimizer in the read

			if(min_strand == 0){
				current_signature.insert(seq + min_loc - kmer_len + 1);
				bin_no = s_mapper->get_bin_id(current_signature.get());
				bins[bin_no]->PutExtendedKmer(seq + min_loc - kmer_len + 1, kmer_len);
			}
			else{
				for(uint32_t rev_i = 0; rev_i < signature_len; rev_i++){
					rev[rev_i] = 3 - (int)seq[min_loc - rev_i];
				}
				
				current_signature.insert(rev);
				bin_no = s_mapper->get_bin_id(current_signature.get());
				bins[bin_no]->PutExtendedKmer(seq + min_loc - kmer_len + 1, kmer_len);
			}

		}
	}

	pmm_reads->free(seq);
	pmm_reads->free(rev);

	return true;
}

//----------------------------------------------------------------------------------
// Process the reads from the given FASTQ file part in small k optimization mode
template<typename COUNTER_TYPE>
bool CSplitter::ProcessReadsSmallK(uchar *_part, uint64 _part_size, ReadType read_type, CSmallKBuf<COUNTER_TYPE>& small_k_buf)
{
	part = _part;
	part_size = _part_size;
	part_pos = 0;

	char *seq;
	uint32 seq_size;
	int omit_next_n_kmers;
	CKmer<1> kmer_str, kmer_rev, kmer_can;
	uint32 i;
	CKmer<1> kmer_mask;
	pmm_reads->reserve(seq);
	kmer_mask.set_n_1(2 * kmer_len);

	uint32 kmer_len_shift = (kmer_len - 1) * 2;

	if (both_strands)
		while (GetSeq(seq, seq_size, read_type))
		{
			if (homopolymer_compressed)
				HomopolymerCompressSeq(seq, seq_size);
			//if (file_type != multiline_fasta)
			//	n_reads++;

			// Init k-mer
			kmer_str.clear();
			kmer_rev.clear();

			// Process first k-1 symbols of a read
			uint32 str_pos = kmer_len_shift - 2;
			uint32 rev_pos = 2;

			omit_next_n_kmers = 0;

			for (i = 0; i < kmer_len - 1; ++i, str_pos -= 2, rev_pos += 2)
			{
				if (seq[i] < 0)
				{
					seq[i] = 0;
					omit_next_n_kmers = i + 1;
				}
				kmer_str.set_2bits(seq[i], str_pos);
				kmer_rev.set_2bits(3 - seq[i], rev_pos);
			}

			// Process next part of a read
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)		// N in a read
				{
					seq[i] = 0;
					omit_next_n_kmers = kmer_len;		// Mark how many symbols to ommit to get the next kmer without any N
				}
				kmer_str.SHL_insert_2bits(seq[i]);
				kmer_str.mask(kmer_mask);
				kmer_rev.SHR_insert_2bits(3 - seq[i], kmer_len_shift);

				// If necessary ommit next symbols
				if (omit_next_n_kmers > 0)
				{
					omit_next_n_kmers--;
					continue;
				}

				// Find canonical kmer representation
				kmer_can = (kmer_str < kmer_rev) ? kmer_str : kmer_rev;

				++small_k_buf.buf[kmer_can.data];
				++total_kmers;
			}
		}
	else
		while (GetSeq(seq, seq_size, read_type))
		{
			if (homopolymer_compressed)
				HomopolymerCompressSeq(seq, seq_size);
			//if (file_type != multiline_fasta)
			//	n_reads++;

			// Init k-mer
			kmer_str.clear();

			// Process first k-1 symbols of a read
			uint32 str_pos = kmer_len_shift - 2;

			omit_next_n_kmers = 0;

			for (i = 0; i < kmer_len - 1; ++i, str_pos -= 2)
			{
				if (seq[i] < 0)
				{
					seq[i] = 0;
					omit_next_n_kmers = i + 1;
				}
				kmer_str.set_2bits(seq[i], str_pos);
			}

			// Process next part of a read
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)		// N in a read
				{
					seq[i] = 0;
					omit_next_n_kmers = kmer_len;		// Mark how many symbols to ommit to get the next kmer without any N
				}
				kmer_str.SHL_insert_2bits(seq[i]);
				kmer_str.mask(kmer_mask);

				// If necessary ommit next symbols
				if (omit_next_n_kmers > 0)
				{
					omit_next_n_kmers--;
					continue;
				}

				++small_k_buf.buf[kmer_str.data];
				++total_kmers;
			}
		}

	pmm_reads->free(seq);
	return true;
}

//----------------------------------------------------------------------------------
// Finish the processing of input file
void CSplitter::Complete()
{
	for (auto& bin : bins)
		if (bin)
			bin->Flush();
}

//************************************************************************************************************
// CWSplitter class - wrapper for multithreading purposes
//************************************************************************************************************
//----------------------------------------------------------------------------------
// Constructor
CWSplitter::CWSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	pq = Queues.part_queue.get();
	bpq = Queues.bpq.get();
	pmm_fastq = Queues.pmm_fastq.get();
	spl = std::make_unique<CSplitter>(Params, Queues);
	spl->InitBins(Params, Queues);
}

//----------------------------------------------------------------------------------
// Execution
void CWSplitter::operator()()
{
	// Splitting parts
	while (!pq->completed())
	{
		uchar *part;
		uint64 size;
		ReadType read_type;
		
		if (pq->pop(part, size, read_type))
		{			
			spl->ProcessReads(part, size, read_type);
			pmm_fastq->free(part);
		}
	}
	spl->Complete();
	bpq->mark_completed();

	spl->GetTotal(n_reads);

	spl.reset();
}

//----------------------------------------------------------------------------------
// Destructor
CWSplitter::~CWSplitter()
{
}

//----------------------------------------------------------------------------------
// Return statistics
void CWSplitter::GetTotal(uint64 &_n_reads)
{
	if (spl)
		spl->GetTotal(n_reads);

	_n_reads = n_reads;
}


//************************************************************************************************************
// CWStatsSplitter class - wrapper for multithreading purposes
//************************************************************************************************************


//----------------------------------------------------------------------------------
// Constructor
CWStatsSplitter::CWStatsSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	spq = Queues.stats_part_queue.get();
	pmm_fastq = Queues.pmm_fastq.get();
	pmm_stats = Queues.pmm_stats.get();
	progressObserver = Params.progressObserver;
	spl = std::make_unique<CSplitter>(Params, Queues);

	signature_len = Params.signature_len;
	pmm_stats->reserve(stats);
	fill_n(stats, (1 << signature_len * 2) + 1, 0);
}

//----------------------------------------------------------------------------------
// Destructor
CWStatsSplitter::~CWStatsSplitter()
{
	pmm_stats->free(stats);
}

//----------------------------------------------------------------------------------
// Execution
void CWStatsSplitter::operator()()
{
	// Splitting parts
	while (!spq->completed())
	{
		uchar *part;
		uint64 size;
		
		ReadType read_type;
		if (spq->pop(part, size, read_type))
		{
			spl->CalcStats(part, size, read_type, stats);
			progressObserver->Step();
			pmm_fastq->free(part);
		}
	}

	spl.reset();
}

//----------------------------------------------------------------------------------
void CWStatsSplitter::GetStats(uint32* _stats)
{
	uint32 size = (1 << signature_len * 2) + 1;
	for (uint32 i = 0; i < size; ++i)
		_stats[i] += stats[i];
}


//************************************************************************************************************
// CWSmallKSplitter class - wrapper for multithreading purposes
//************************************************************************************************************


//----------------------------------------------------------------------------------
// Constructor
template <typename COUNTER_TYPE> CWSmallKSplitter<COUNTER_TYPE>::CWSmallKSplitter(CKMCParams &Params, CKMCQueues &Queues)
{
	pq = Queues.part_queue.get();
	pmm_fastq = Queues.pmm_fastq.get();
	pmm_small_k = Queues.pmm_small_k_buf.get();
	kmer_len = Params.kmer_len;
	spl = std::make_unique<CSplitter>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
template <typename COUNTER_TYPE> CWSmallKSplitter<COUNTER_TYPE>::~CWSmallKSplitter()
{
}

//----------------------------------------------------------------------------------
// Execution
template <typename COUNTER_TYPE> void CWSmallKSplitter<COUNTER_TYPE>::operator()()
{
	pmm_small_k->reserve(small_k_buf.buf);
	memset(small_k_buf.buf, 0, (1ull << 2 * kmer_len) * sizeof(*small_k_buf.buf));

	// Splitting parts
	while (!pq->completed())
	{
		uchar *part;
		uint64 size;
		ReadType read_type;
		if (pq->pop(part, size, read_type))
		{
			spl->ProcessReadsSmallK(part, size, read_type, small_k_buf);
			pmm_fastq->free(part);
		}
	}
	spl->Complete();

	spl->GetTotal(n_reads);
	total_kmers = spl->GetTotalKmers();
	spl.reset();
}

//----------------------------------------------------------------------------------
// Return statistics
template <typename COUNTER_TYPE> void CWSmallKSplitter<COUNTER_TYPE>::GetTotal(uint64 &_n_reads)
{
	if (spl)
		spl->GetTotal(n_reads);
	_n_reads = n_reads;
}




//************************************************************************************************************
// CWSplitter class - wrapper for multithreading purposes
//************************************************************************************************************
//----------------------------------------------------------------------------------
// Constructor
CWEstimateOnlySplitter::CWEstimateOnlySplitter(CKMCParams& Params, CKMCQueues& Queues)
{
	pq = Queues.part_queue.get();
	pmm_fastq = Queues.pmm_fastq.get();
	spl = std::make_unique<CSplitter>(Params, Queues);
}

//----------------------------------------------------------------------------------
// Execution
void CWEstimateOnlySplitter::operator()()
{
	// Splitting parts
	while (!pq->completed())
	{
		uchar* part;
		uint64 size;
		ReadType read_type;
		
		if (pq->pop(part, size, read_type))
		{
			spl->ProcessReadsOnlyEstimate(part, size, read_type);
			pmm_fastq->free(part);
		}
	}

	spl->GetTotal(n_reads);

	spl.reset();
}

//----------------------------------------------------------------------------------
// Destructor
CWEstimateOnlySplitter::~CWEstimateOnlySplitter()
{
}

//----------------------------------------------------------------------------------
// Return statistics
void CWEstimateOnlySplitter::GetTotal(uint64& _n_reads)
{
	if (spl)
		spl->GetTotal(n_reads);

	_n_reads = n_reads;
}





//instantiate some templates
template bool CSplitter::ProcessReadsSmallK(uchar *_part, uint64 _part_size, ReadType read_type, CSmallKBuf<uint32>& small_k_buf);
template bool CSplitter::ProcessReadsSmallK(uchar *_part, uint64 _part_size, ReadType read_type, CSmallKBuf<uint64>& small_k_buf);
template class CWSmallKSplitter<uint32>;
template class CWSmallKSplitter<uint64>;

// ***** EOF
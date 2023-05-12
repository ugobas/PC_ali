struct rank{
  float score;
  int index;
  struct rank *next;
};
int Rank(float score, struct rank *rank, 
	 int *ranked, int r_max,
	 struct rank **First_rank, struct rank **Last_rank);

#include "rank.h"
#include <stddef.h>

int Rank(float score, struct rank *rank, 
	 int *ranked, int r_max,
	 struct rank **First_rank, struct rank **Last_rank)
{
  int k;
  if(*ranked < r_max){k=*ranked; (*ranked)++;}
  else{
    if(score < (*Last_rank)->score)return(-1);
    k=(*Last_rank)->index;
  }
  struct rank *toberanked=rank+k;
  toberanked->index=k;
  toberanked->score=score;
  if(*ranked==1){
    *First_rank=toberanked;
    *Last_rank=toberanked;
    toberanked->next=NULL;
    return(k);
  }
  struct rank *tmp=*First_rank, *prev=NULL, *next;
  while(tmp && tmp->score > toberanked->score){
    prev=tmp; tmp=tmp->next;
  }
  if(prev==NULL){
    next=*First_rank;
    *First_rank=toberanked;
  }else{
    next=prev->next;
    prev->next=toberanked;
  }
  toberanked->next=next;
  if(next==NULL){*Last_rank=toberanked;}
  return(k);
}

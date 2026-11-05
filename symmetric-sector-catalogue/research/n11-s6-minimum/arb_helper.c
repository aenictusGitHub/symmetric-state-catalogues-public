/* Arithmetic-only accelerator. Opaque FLINT objects are accessed exclusively
   through exported FLINT functions, without assuming their internal layout. */
#include <stdlib.h>
#include <stdint.h>
extern void *_acb_vec_init(long);
extern void _acb_vec_clear(void *,long);
extern void *_acb_vec_entry_ptr(void *,long);
extern void acb_mat_init(void *,long,long);
extern void acb_mat_clear(void *);
extern void *acb_mat_entry_ptr(void *,long,long);
extern void acb_mat_mul(void *,const void *,const void *,long);
extern void *acb_real_ptr(void *);
extern void *acb_imag_ptr(void *);
extern void acb_zero(void *);
extern void acb_set(void *,const void *);
extern void acb_set_si_si(void *,long,long);
extern void acb_mul_2exp_si(void *,const void *,long);
extern void acb_add(void *,const void *,const void *,long);
extern void acb_mul(void *,const void *,const void *,long);
extern void acb_addmul(void *,const void *,const void *,long);
extern void acb_mul_si(void *,const void *,long,long);
extern void acb_conj(void *,const void *);
extern int arb_set_str(void *,const char *,long);
extern char *arb_get_str(const void *,long,unsigned long);
extern void flint_free(void *);
static void *ve(void*v,long i){return _acb_vec_entry_ptr(v,i);}
void *newmat(long r,long c){void *m=calloc(1,256);acb_mat_init(m,r,c);return m;}
void delmat(void*m){acb_mat_clear(m);free(m);}
void *newvec(long n){return _acb_vec_init(n);}
void delvec(void*v,long n){_acb_vec_clear(v,n);}
void *matentry(void*m,long r,long c){return acb_mat_entry_ptr(m,r,c);}
void *vecentry(void*v,long i){return ve(v,i);}
int setvalue(void*z,const char*re,const char*im,long prec){return arb_set_str(acb_real_ptr(z),re,prec)||arb_set_str(acb_imag_ptr(z),im,prec);}
char *realvalue(void*z,long digits){return arb_get_str(acb_real_ptr(z),digits,0);}
void freestring(void*s){flint_free(s);}
void eval_source(void*out,void*mix,const int64_t*digits,long d,long ncol,long prec){
 void *raw=newmat(d,ncol),*temp=newvec(2);
 for(long i=0;i<d;i++)for(long j=0;j<ncol;j++){
  long k=6*(i*ncol+j);void*z=matentry(raw,i,j);
  acb_set_si_si(z,digits[k+2],digits[k+5]);acb_mul_2exp_si(z,z,32);
  acb_set_si_si(ve(temp,0),digits[k+1],digits[k+4]);acb_add(z,z,ve(temp,0),prec);acb_mul_2exp_si(z,z,32);
  acb_set_si_si(ve(temp,0),digits[k],digits[k+3]);acb_add(z,z,ve(temp,0),prec);
 }
 acb_mat_mul(out,mix,raw,prec);delmat(raw);delvec(temp,2);
}
/* ops = ci, source-block, source-channel, source-component, fundamental-index,
   output-component. Source dimensions are supplied as p,q; output p,q is fixed. */
void eval_block(void*out,void**sources,const int64_t*sourcepq,const int64_t*ops,
 void*coeff,long nops,void*mix,void*logical,void*weights,const int64_t*point,
 long d,long L,long p,long q,const int64_t*piv,long npiv,void*answer,long prec){
 long np=(p+1)*(q+1);void *raw=newmat(d,(L+1)*np),*F=newmat(d,(L+1)*np),*E=newmat(d,(L+1)*np),*tmp=newvec(5),*xy=newvec(24);
 for(long b=0;b<12;b++){
  acb_set_si_si(ve(xy,b),point[2*b],point[2*b+1]);long sg=b%2?-1:1;
  acb_set_si_si(ve(xy,12+b),sg*point[2*(11-b)],-sg*point[2*(11-b)+1]);
 }
 for(long h=0;h<nops;h++){
  long ci=ops[6*h],ib=ops[6*h+1],j=ops[6*h+2],a=ops[6*h+3],b=ops[6*h+4],r=ops[6*h+5];
  long ps=sourcepq[2*ib],qs=sourcepq[2*ib+1],nsp=(ps+1)*(qs+1);int f=(ps+1==p);
  for(long i=0;i<=ps;i++)for(long t=0;t<=qs;t++){
   acb_mul(ve(tmp,0),matentry(sources[ib],j,a*nsp+i*(qs+1)+t),ve(coeff,h),prec);
   if(f){
    acb_addmul(matentry(raw,ci,r*np+i*(q+1)+t),ve(tmp,0),ve(xy,b),prec);
    acb_addmul(matentry(raw,ci,r*np+(i+1)*(q+1)+t),ve(tmp,0),ve(xy,12+b),prec);
   }else{
    acb_addmul(matentry(raw,ci,r*np+i*(q+1)+t),ve(tmp,0),ve(xy,12+b),prec);
    acb_mul_si(ve(tmp,1),ve(tmp,0),-1,prec);
    acb_addmul(matentry(raw,ci,r*np+i*(q+1)+t+1),ve(tmp,1),ve(xy,b),prec);
   }
  }
 }
 acb_mat_mul(F,mix,raw,prec);
 for(long a=0;a<d;a++)for(long r=0;r<=L;r++)for(long t=0;t<np;t++)for(long u=0;u<np;u++)
  acb_addmul(matentry(E,a,r*np+t),matentry(F,a,r*np+u),matentry(logical,u,t),prec);
 /* Identity Gram contribution and the selected correction-entry features. */
 for(long a=0;a<d;a++)for(long r=0;r<=L;r++)for(long t=0;t<np;t++){
  void*z=matentry(E,a,r*np+t);acb_conj(ve(tmp,0),z);acb_mul(ve(tmp,1),z,ve(tmp,0),prec);
  acb_addmul(ve(answer,0),ve(tmp,1),ve(weights,r*np+t),prec);
 }
 for(long h=0;h<npiv;h++){
  long index=piv[3*h],a=piv[3*h+1],b=piv[3*h+2];void*res=ve(answer,1+index);
  for(long r=0;r<=L;r++)for(long t=0;t<np;t++){
   acb_conj(ve(tmp,0),matentry(E,b,r*np+t));acb_mul(ve(tmp,1),matentry(E,a,r*np+t),ve(tmp,0),prec);
   if(a!=b)acb_mul_si(ve(tmp,1),ve(tmp,1),2,prec);
   acb_addmul(res,ve(tmp,1),ve(weights,r*np+t),prec);
  }
 }
 delmat(raw);delmat(F);delmat(E);delvec(tmp,5);delvec(xy,24);
}

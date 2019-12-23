#include<stdio.h>
#include<gmp.h>
#include<stdlib.h>
#include<math.h>

mpz_t p,a,b;

void print_point(mpz_t x,mpz_t y)
{
  	gmp_printf("x3 : %Zd ",x);
  	gmp_printf("y3 : %Zd\n",y);
}

mpz_t * scan_point()
{
	mpz_t *Coordinates=malloc(sizeof(mpz_t)*2);
	mpz_init(Coordinates[0]);
	mpz_init(Coordinates[1]);
	printf("Enter x and y\n");
	gmp_scanf("%Zd",Coordinates[0]);
	gmp_scanf("%Zd",Coordinates[1]);

	return Coordinates;
}


mpz_t * add(mpz_t x1,mpz_t y1,mpz_t x2,mpz_t y2)
{
	//sum between two points
	mpz_t *Points;
	Points=malloc(sizeof(mpz_t)*2);
	mpz_init(Points[0]);
	mpz_init(Points[1]);

	if(mpz_cmp_si(x1,-1)==0 && mpz_cmp_si(y1,-1)==0)
	{
		mpz_set(Points[0],x2);
		mpz_set(Points[1],y2);
	}

	else if(mpz_cmp_si(x2,-1)==0 && mpz_cmp_si(y2,-1)==0)
	{
		mpz_set(Points[0],x1);
		mpz_set(Points[1],y1);
	}

	//point doubling
	else if(mpz_cmp(x1,x2)==0 && mpz_cmp(y1,y2)==0)
	{
		//initialize
		mpz_t x3,y3,lambda,temp,flag;
		mpz_init(x3);
		mpz_init(y3);
		mpz_init(temp);
		mpz_init(flag);
		mpz_init(lambda);

		//calculating lambda
		mpz_mul_ui(temp,y1,2);
		mpz_invert(temp,temp,p);
		mpz_mul(flag,x2,x2);
		mpz_mul_ui(flag,flag,3);
		mpz_add(flag,flag,a);
		mpz_mul(lambda,flag,temp);
		mpz_mod(lambda,lambda,p);

		//calculating x3
		mpz_mul(flag,lambda,lambda);
		mpz_mul_ui(temp,x1,2);
		mpz_sub(x3,flag,temp);
		mpz_mod(x3,x3,p);
		mpz_set(Points[0],x3);

		//calculating y3
		mpz_sub(temp,x2,x3);
		mpz_mul(flag,temp,lambda);
		mpz_sub(y3,flag,y1);
		mpz_mod(y3,y3,p);
		mpz_set(Points[1],y3);


		mpz_clear(x3);
		mpz_clear(y3);
		mpz_clear(lambda);
		mpz_clear(temp);
		mpz_clear(flag);

	}

	//adding two points with diff x coordinates
	else if(mpz_cmp(x1,x2)!=0)
	{
		mpz_t x3,y3,lambda,temp,flag;
		mpz_init(x3);
		mpz_init(y3);
		mpz_init(temp);
		mpz_init(flag);
		mpz_init(lambda);
		

		//calculating lambda
		mpz_sub(temp,y2,y1);
		mpz_sub(flag,x2,x1);
		mpz_invert(flag,flag,p);
		mpz_mul(lambda,temp,flag);
		mpz_mod(lambda,lambda,p);

		//calculating x3
		mpz_mul(temp,lambda,lambda);
		mpz_sub(flag,temp,x1);
		mpz_sub(x3,flag,x2);
		mpz_mod(x3,x3,p);
		mpz_set(Points[0],x3);

		//calculating y3
		mpz_sub(temp,x1,x3);
		mpz_mul(flag,lambda,temp);
		mpz_sub(y3,flag,y1);
		mpz_mod(y3,y3,p);
		mpz_set(Points[1],y3);

	
		mpz_clear(x3);
		mpz_clear(y3);
		mpz_clear(lambda);
		mpz_clear(temp);
		mpz_clear(flag);


	}
	//same x and different y, or y==0
	else if(mpz_cmp(x1,x2)==0 && ((mpz_cmp(y1,y2)!=0) || (mpz_cmp_ui(y1,0)==0 && mpz_cmp_ui(y2,0)==0)))
	{
		mpz_set_si(Points[0],-1);
		mpz_set_si(Points[1],-1);	
	}

return Points;
	
}

mpz_t * negation(mpz_t x,mpz_t y)
{

	mpz_t n;
	mpz_init(n);
	mpz_t *Points;
	Points=malloc(sizeof(mpz_t)*2);
	mpz_init(Points[0]);
	mpz_init(Points[1]);
	mpz_set(Points[0],x);
	if(mpz_cmp_si(y,-1)==0 && mpz_cmp_si(x,-1)==0)
		mpz_set(Points[1],y);
	else
	{
		mpz_neg(n,y);
		mpz_mod(n,n,p);
		mpz_set(Points[1],n);
	}

	
	mpz_clear(n);
	return Points;

}

mpz_t * subtract(mpz_t x,mpz_t y,mpz_t m,mpz_t n)
{
	mpz_t *pt=negation(m,n);
  	mpz_t *point=add(x,y,pt[0],pt[1]);
	mpz_clear(pt[0]);
  	mpz_clear(pt[1]);
  	free(pt);
  	return point;
}

mpz_t* shanks_algo(mpz_t A)
{
	mpz_t *Y=malloc(sizeof(mpz_t));
	mpz_t s,z,x,temp,temp3,two,q;
	mpz_t t,u,k,B,y,m;
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	mpz_init(s);
	mpz_init(t);
	mpz_init(z);
	mpz_init(u);
	mpz_init(k);
	mpz_init(x);
	mpz_init(B);
	mpz_init(y);
	mpz_init(m);
	mpz_init(q);
	mpz_init(two);
	mpz_init(temp);
	mpz_init(temp3);
	mpz_init(Y[0]);
	mpz_set_ui(s,0);
	mpz_sub_ui(q,p,1);
	mpz_set(t,q);
	mpz_set_ui(two,2);

	while(mpz_even_p(t)!=0)
	{
		mpz_divexact_ui(t,t,2);
		mpz_add_ui(s,s,1);
	}

	do
	{
		mpz_urandomm(u,state,p);

	}while(mpz_legendre(u,p)!=-1);


	mpz_set(k,s);	//k=s
	mpz_powm(z,u,t,p); //z=u^t mod p

	mpz_add_ui(temp,t,1);
	mpz_divexact_ui(temp,temp,2); //calc (t+1)/2

	mpz_powm(x,A,temp,p);  // x= a^((t+1)/2) mod p
	mpz_powm(B,A,t,p); //b=a^t mod p

	if(mpz_cmp_ui(B,0)==0)
	{
		mpz_set(Y[0],B);
	}
	
  else
  {

	while(mpz_cmp_ui(B,1)!=0)
	{
		mpz_set_ui(m,1);
		mpz_set_ui(temp,1);
		while(mpz_cmp(m,k)<0)
		{
			mpz_mul_ui(temp,temp,2);
			mpz_powm(temp3,B,temp,p);
			if(mpz_cmp_ui(temp3,1)==0)
			{
				break;
			}
			else
			{
				mpz_add_ui(m,m,1);
			}
		}

		mpz_sub(temp3,k,m);
		mpz_sub_ui(temp3,temp3,1);
		mpz_sub_ui(q,p,1);
		mpz_powm(temp,two,temp3,q); //2^k-m-1
		mpz_powm(y,z,temp,p); //y=z^2^k-m-1 mod p
		mpz_powm_ui(z,y,2,p); //z=y^2 mod p
		mpz_mul(B,B,z);//B=Bz mod p
		mpz_mod(B,B,p);
		mpz_mul(x,x,y); //x=xy mod p
		mpz_mod(x,x,p);
		mpz_set(k,m);
		
	}

	mpz_set(Y[0],x);
  }

	mpz_clear(s);
	mpz_clear(t);
	mpz_clear(u);
	mpz_clear(z);
	mpz_clear(k);
	mpz_clear(y);
	mpz_clear(m);
	mpz_clear(temp);
	mpz_clear(temp3);
	gmp_randclear(state);
	mpz_clear(x);
	mpz_clear(two);
	mpz_clear(B);
	mpz_clear(q);

	return Y;

		
}

mpz_t * generate(mpz_t x)
{
	mpz_t temp,ax;
	mpz_t *val=malloc(sizeof(mpz_t));
	mpz_init(ax);
	mpz_init(temp);
	mpz_init(val[0]);
	mpz_set_si(val[0],-1);
	mpz_pow_ui(temp,x,3);
	mpz_mul(ax,a,x);
	mpz_add(temp,temp,ax);
	mpz_add(temp,temp,b);
	mpz_mod(temp,temp,p);
	mpz_clear(ax);
	if(mpz_legendre(temp,p)==-1)
		return val;
	else
		val=shanks_algo(temp);
	
	return val;

}

mpz_t * point_multiplication(mpz_t x,mpz_t y,mpz_t k)
{
  	unsigned long int *bits=(unsigned long int*)(malloc(sizeof(int)*512));
  	int c=0;
  	while(mpz_cmp_ui(k,0)!=0)
  	{
  		bits[c]=mpz_fdiv_q_ui(k,k,2);
  		c++;
  	}

  	mpz_t *P=malloc(sizeof(mpz_t)*2);
  	mpz_init(P[0]);
  	mpz_init(P[1]);
  	mpz_set_si(P[0],-1);
  	mpz_set_si(P[1],-1);
  	for(int j=c-1;j>=0;j--)
  	{
  		P=add(P[0],P[1],P[0],P[1]);
  		if(bits[j]==1)
  			P=add(P[0],P[1],x,y);
  	}


  	free(bits);
  	return P;
  				
}

void read_from_file(FILE *FP)
{
	mpz_inp_str(p,FP,10); 
	mpz_inp_str(a,FP,10);
	mpz_inp_str(b,FP,10);
}


int main(int argc,char *argv[])
{
	FILE *fp=fopen("Ecurve.txt","r");
	mpz_t x,y,m,n;
	mpz_init(p);
	mpz_init(a);
	mpz_init(b);
	mpz_init(x);
	mpz_init(y);
	mpz_init(m);
	mpz_init(n);
	
	//read domain parameters from file:
	read_from_file(fp);
	int choice=0;
	printf("1 for Adding two points\n" 
		"2 for Negation of a point\n"
		"3 for Subtraction of 2 Points\n"
		"4 for Point Multiplication\n"
		"5 to find y for given x\n");

	scanf("%d",&choice);

		if(choice==1)
		{
			mpz_t *C1=scan_point();
  			mpz_t *C2=scan_point();
  			mpz_t *P=add(C1[0],C1[1],C2[0],C2[1]);
  			print_point(P[0],P[1]);
  			free(P);
  			free(C1);
  			free(C2);
  		}


  		else if(choice==2)
  		{
  			mpz_t *c=scan_point();
  			mpz_t *Pt=negation(c[0],c[1]);
  			print_point(Pt[0],Pt[1]);
  			free(Pt);
  			free(c);	
  		}
				
  		else if(choice==3)	
  		{
  			mpz_t *c1=scan_point();
  			mpz_t *c2=scan_point();
  			mpz_t *point=subtract(c1[0],c1[1],c2[0],c2[1]);
  			print_point(point[0],point[1]);
  			free(point);
  			free(c1);
  			free(c2);
  			

  		}

  		else if(choice==4)
  		{
  			mpz_t *C=scan_point();
  			printf("Enter value of k\n");
  			mpz_t k;
  			mpz_init(k);
  			gmp_scanf("%Zd",k);
  			mpz_t *PM=point_multiplication(C[0],C[1],k);
  			print_point(PM[0],PM[1]);
  			mpz_clear(k);
  			free(PM);
  			free(C);

  		}		

  		else if(choice == 5)
  		{
  			mpz_t *val=malloc(sizeof(mpz_t));
  			mpz_init(val[0]);
  			printf("Enter the value of x to find its y\n");
  			gmp_scanf("%Zd",x);
  			val=generate(x);
  			if(mpz_cmp_si(val[0],-1)==0)
  					printf("Points dont exist on curve\n");
  			else
  			{
  				gmp_printf("New y: %Zd\n",val[0]);
  			}
  			mpz_clear(val[0]);
  			free(val);
  		} 	
  				

  		else
  		{
  			printf("Not available! Try again\n");
  		}		
  	  

	mpz_clear(p);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(m);
	mpz_clear(n);
	fclose(fp);
}
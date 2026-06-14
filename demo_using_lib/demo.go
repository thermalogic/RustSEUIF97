// The Go example to call the shared library: cdel
// make -f mgo.mk
package main
/*
double pt(double p, double t, int o_id);
double pt2s(double p, double t);
*/
import "C"
import "fmt"

func main() {
	p:=16.14
	t:=512.4
	oid:=4
	h := C.pt(C.double(p),C.double(t),C.int(oid))
	s:= C.pt2s(C.double(p),C.double(t))
	fmt.Printf("(p,t)=(%.2f, %.2f), h=%.2f s=%.4f\n",p,t,float64(h),float64(s))
}
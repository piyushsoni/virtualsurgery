class SymbolEncoderStream extends FilterOutputStream {
  int[] C = {0};
  int[] R = {1, 0};
  int[] S = {1, 1, 0};
  int[] L = {1, 1, 1, 0};
  int[] M = {1, 1, 1, 1, 0};
  int[] E = {1, 1, 1, 1, 1};
  
  int[] buf = new int[8];
  int bptr = 0;
  
  SymbolEncoderStream(OutputStream os) { super(os); }
  
  int intVal(int[] b) {return b[0]*128 + b[1]*64 + b[2]*32 + b[3]*16 + b[4]*8 + b[5]*4 + b[6]*2 + b[7]; }
  
  void doWrite(int[] b) throws IOException { // add symbol code to buffer, write buffer contents to underlying stream if full
     for (int i = 0; i < b.length; i++) {
         buf[bptr++] = b[i];
         if (bptr==8) {out.write(intVal(buf)); bptr = 0;}
     }
  }
  
  void flush() throws IOException {    // append buffer content with zeros and write it to underlying stream
     for (int i = bptr; i < 8; i++) buf[i] = 0; 
     out.write(intVal(buf)); bptr = 0;
     out.flush();
  }
  
  void write(int b) throws IOException {
     switch (b) {
        case 'C': doWrite(C); break;
        case 'L': doWrite(L); break;
        case 'R': doWrite(R); break;
        case 'S': doWrite(S); break;
        case 'M': doWrite(M); break;
        case 'E': doWrite(E); break;
     } 
  }
  
}

class SymbolDecoderStream extends FilterInputStream {
  
    int[] buf = new int[16];
    int bptr = 0;
  
    SymbolDecoderStream(InputStream is) {super(is); }
    
    boolean doRead() throws IOException { // read next byte into buffer
       int b = in.read();  
       if (b == -1) return false;  
       for (int i = 7; i >= 0; i--) { buf[bptr+i] = b % 2; b /= 2; }
       bptr += 8;
       return true;
    }
    
    int parseSymbol() {    // try to parse a symbol from the current buffer contents
       int s = -1; int l = 0;
       if (buf[0] == 0 && bptr > 0) { s='C'; l=1; }
       else if (buf[1] == 0 && bptr > 1) {s='R'; l=2;}
       else if (buf[2] == 0 && bptr > 2) {s='S'; l=3;}
       else if (buf[3] == 0 && bptr > 3) {s='L'; l=4;}
       else if (buf[4] == 0 && bptr > 4) {s='M'; l=5;}
       else if (buf[4] == 1 && bptr > 4) {s='E'; l=5;}
       for (int i = 0; i < bptr - l; i++) buf[i] = buf[i+l];
       bptr -= l;
       return s; 
    }
    
    int read() throws IOException {
        int s;
        while (-1 == (s = parseSymbol())) {if (!doRead()) return -1;}  // read bytes as long as we cannot decode any symbols
        return s;       
    }
}



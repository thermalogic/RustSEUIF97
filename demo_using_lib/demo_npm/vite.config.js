import { defineConfig } from 'vite';

export default defineConfig({
  server: {
    port: 3000
  },
  assetsInclude: ['**/*.wasm'],
  optimizeDeps: {
    exclude: ['seuif97']
  }
});
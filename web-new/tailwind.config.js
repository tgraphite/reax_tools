/** @type {import('tailwindcss').Config} */
module.exports = {
  content: [
    './*.html',
    './src/**/*.{js,html}',
  ],
  darkMode: 'class',
  theme: {
    extend: {
      colors: {
        primary: {
          50: '#ebf5fb',
          100: '#d6eaf8',
          200: '#aed6f1',
          300: '#85c1e9',
          400: '#5dade2',
          500: '#3498db',
          600: '#2980b9',
          700: '#2471a3',
          800: '#1f618d',
          900: '#1a5276',
        },
        secondary: {
          50: '#f4f6f7',
          100: '#e5e8e8',
          200: '#ccd1d1',
          300: '#b2babb',
          400: '#99a3a4',
          500: '#7f8c8d',
          600: '#616a6b',
          700: '#515a5a',
          800: '#424949',
          900: '#2c3e50',
        },
        success: {
          500: '#27ae60',
          600: '#229954',
        },
        warning: {
          500: '#f39c12',
          600: '#d68910',
        },
        danger: {
          500: '#e74c3c',
          600: '#c0392b',
        },
      },
      fontFamily: {
        sans: ['Inter', 'system-ui', '-apple-system', 'BlinkMacSystemFont', 'Segoe UI', 'Roboto', 'sans-serif'],
        mono: ['JetBrains Mono', 'Fira Code', 'Consolas', 'Monaco', 'monospace'],
      },
      animation: {
        'pulse-slow': 'pulse 3s cubic-bezier(0.4, 0, 0.6, 1) infinite',
        'fade-in': 'fadeIn 0.5s ease-out',
        'slide-up': 'slideUp 0.3s ease-out',
      },
      keyframes: {
        fadeIn: {
          '0%': { opacity: '0' },
          '100%': { opacity: '1' },
        },
        slideUp: {
          '0%': { transform: 'translateY(10px)', opacity: '0' },
          '100%': { transform: 'translateY(0)', opacity: '1' },
        },
      },
    },
  },
  plugins: [
    require('@tailwindcss/forms'),
  ],
}

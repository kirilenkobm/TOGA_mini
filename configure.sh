#!/usr/bin/env bash
set -e

echo "🔧 TOGA-mini Configuration Script"
echo "=================================="

# Check if we're in a conda environment
if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "⚠️  Warning: Not in a conda environment. Consider running:"
    echo "   conda env create -f environment.yaml"
    echo "   conda activate TOGA-mini"
    echo ""
fi

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check required tools
echo "🔍 Checking dependencies..."
missing_deps=()

if ! command_exists cmake; then
    missing_deps+=("cmake")
fi

if ! command_exists make; then
    missing_deps+=("make")
fi

if ! command_exists gcc; then
    missing_deps+=("gcc")
fi

if [ ${#missing_deps[@]} -ne 0 ]; then
    echo "❌ Missing dependencies: ${missing_deps[*]}"
    echo "Please install them with: conda install ${missing_deps[*]}"
    exit 1
fi

echo "✅ All dependencies found"

# Build C modules
echo ""
echo "🏗️  Building C modules..."
./build_c.sh

# Train models
echo ""
echo "🧠 Training classification models..."
if [ -f "chain_class_models/train_toga_chain_class_model.py" ]; then
    cd chain_class_models
    python train_toga_chain_class_model.py
    cd ..
    echo "✅ Models trained successfully"
else
    echo "⚠️  Warning: Model training script not found, skipping..."
fi

# Verify installation
echo ""
echo "🔍 Verifying installation..."
python -c "
import sys
import os
sys.path.insert(0, '.')

try:
    from toga_modules.common import to_log
    print('✅ Python modules import successfully')
except ImportError as e:
    print(f'❌ Python module import failed: {e}')
    sys.exit(1)

# Check C libraries
util_c_lib = 'util_c/lib'
if os.path.exists(util_c_lib):
    libs = os.listdir(util_c_lib)
    if libs:
        print(f'✅ C libraries found: {libs}')
    else:
        print('❌ No C libraries found')
        sys.exit(1)
else:
    print('❌ C library directory not found')
    sys.exit(1)

# Check models
models_dir = 'chain_class_models'
required_models = ['se_model.dat', 'me_model.dat']
missing_models = []
for model in required_models:
    if not os.path.exists(os.path.join(models_dir, model)):
        missing_models.append(model)

if missing_models:
    print(f'⚠️  Missing models: {missing_models}')
    print('   You may need to train them manually')
else:
    print('✅ All required models found')
"

echo ""
echo "🎉 Configuration complete!"
echo ""
echo "📖 Usage:"
echo "   ./toga_mini.py --help"
echo ""
echo "📝 For nextflow parallel execution, make sure to:"
echo "   1. Use the same conda environment on all nodes"
echo "   2. Set appropriate nextflow config for your cluster"
echo "   3. Consider using conda-pack for environment distribution" 
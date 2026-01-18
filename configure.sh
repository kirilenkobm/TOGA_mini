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

# Check Python environment
echo "🔍 Checking Python environment..."
python -c "
import sys
print(f'✅ Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}')

try:
    import pyrion
    print('✅ pyrion package available')
except ImportError:
    print('❌ pyrion package not found. Install with: pip install pyrion')
    sys.exit(1)

try:
    import numpy, pandas, joblib, xgboost, sklearn
    print('✅ ML dependencies available')
except ImportError as e:
    print(f'❌ Missing ML dependency: {e}')
    sys.exit(1)
"

# Train models
echo ""
echo "🧠 Training classification models..."
if [ -f "chain_classification_models/train_toga_chain_class_model.py" ]; then
    cd chain_classification_models
    python train_toga_chain_class_model.py
    cd ..
    echo "✅ Models trained successfully"
else
    echo "⚠️  Warning: Model training script not found, skipping..."
fi

echo ""
echo "🎉 Configuration complete!"
echo ""
echo "📖 Usage:"
echo "   python toga_mini.py --help"
echo ""
echo "📝 Quick test:"
echo "   ./quick_test.sh" 
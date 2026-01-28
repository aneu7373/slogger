FROM python:3.11-slim

WORKDIR /app

# Needed for git submodules + building some Python deps if necessary
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
  && rm -rf /var/lib/apt/lists/*

# Copy project first (so Docker cache works better)
COPY pyproject.toml /app/pyproject.toml
COPY src /app/src
COPY tests /app/tests

# If you add codonopt as a submodule at third_party/codonopt, copy it in:
COPY third_party /app/third_party

# Install codonopt (editable) if present
RUN if [ -d /app/third_party/codonopt ]; then \
      pip install --no-cache-dir -e /app/third_party/codonopt; \
    fi

# Install slogger
RUN pip install --no-cache-dir -e /app

ENTRYPOINT ["slogger"]

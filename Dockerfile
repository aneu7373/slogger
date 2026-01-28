FROM python:3.11-slim

WORKDIR /app

# Needed for git submodules (and common build sanity)
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
  && rm -rf /var/lib/apt/lists/*

# Copy project
COPY pyproject.toml /app/pyproject.toml
COPY src /app/src
COPY tests /app/tests

# Copy third_party (codonopt submodule lives here)
COPY third_party /app/third_party

# Make codonopt importable without packaging (no setup.py/pyproject.toml required)
# codonopt package dir is /app/third_party/codonopt/codonopt
ENV PYTHONPATH="/app/third_party/codonopt:${PYTHONPATH}"

# Install codonopt dependencies if present
RUN if [ -f /app/third_party/codonopt/requirements.txt ]; then \
      pip install --no-cache-dir -r /app/third_party/codonopt/requirements.txt; \
    fi

# Install slogger
RUN pip install --no-cache-dir -e /app

ENTRYPOINT ["slogger"]
